import re
import os
import subprocess
from typing import List
from abc import ABCMeta, abstractmethod
import numpy as np
from openbabel import openbabel as ob
from openbabel import pybel


# Copyright (c) 2020 yumn-dfx
# This software is released under the MIT License.
# https://opensource.org/licenses/mit-license.php


class CalcMethods:
    def __init__(self, chkdir='', link='', root='# b3lyp/6-31G*', addsec=''):
        self.chkdir = chkdir
        self.linkop = link
        self.method = root
        self.additional = addsec


class Molecule:
    def __init__(self, pybelmol: pybel.Molecule, name, cm):
        self.pybelobj = pybelmol
        self.name = name
        self.cm = cm
        self.energy = 0
        self.force = np.zeros(len(pybelmol.atoms) * 3)
        self.position = np.array([atom.coords for atom in pybelmol.atoms]).flatten()


class GaussianIO:
    # change the following value to 'g16' if you want to run with Gaussian 16
    gau_cmd = 'g09'

    def __init__(self, method=CalcMethods()):
        self.method = method

    @staticmethod
    def get_energy(outfile):
        energy = -1
        with open(outfile) as fst:
            line = fst.readline()
            while line:
                m = re.search(r'\s*SCF Done.+\s(-?\d+\.\d+)', line)
                if m:
                    energy = float(m.group(1))
                    break
                line = fst.readline()
        return energy
        # To add: check normal termination

    @staticmethod
    def get_force(outfile):
        force_vec = []
        with open(outfile) as fst:
            line = fst.readline()
            while line:
                if 'Forces (Hartrees/Bohr)' in line:
                    fst.readline()
                    fst.readline()
                    line = fst.readline()
                    mm = re.search(r'^\s*?\d+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$', line)
                    while mm:
                        xyz = [float(mm.group(1)), float(mm.group(2)), float(mm.group(3))]
                        force_vec.extend(xyz)
                        line = fst.readline()
                        mm = re.search(r'^\s*?\d+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$', line)
                    break
                line = fst.readline()
        return np.array(force_vec)

    def write_com(self, name, pybelmol, cm, method_additional=''):
        file = name + '.gjf'
        with open(file, mode='w') as fr:
            fr.write('%chk=' + self.method.chkdir + '_' + name + '.chk\n')
            if self.method.linkop != '':
                fr.write(self.method.linkop)
            fr.write(self.method.method + ' ' + method_additional + "\n\n" + "title" + "\n\n")
            fr.write(cm + "\n")
            coord = pybelmol.write("xyz").splitlines()
            coord = coord[2:]
            fr.write('\n'.join(coord) + "\n\n")
            fr.write(self.method.additional + "\n\n")

    def update_energy(self, mol: Molecule, add_calc=''):
        # update only energy
        infile = mol.name + ".gjf"
        outfile = mol.name + ".out"
        run_cmd = f'{self.gau_cmd} < {infile} > {outfile}'
        self.write_com(mol.name, mol.pybelobj, mol.cm, add_calc)
        subprocess.call(run_cmd, shell=True)
        mol.energy = self.get_energy(outfile)

    def update_force(self, mol: Molecule, add_calc=''):
        # update force
        infile = mol.name + ".gjf"
        outfile = mol.name + ".out"
        run_cmd = f'{self.gau_cmd} < {infile} > {outfile}'
        self.write_com(mol.name, mol.pybelobj, mol.cm, add_calc)
        subprocess.call(run_cmd, shell=True)
        mol.energy = self.get_energy(outfile)
        mol.force = self.get_force(outfile)


class Rotate:
    @staticmethod
    def rotate_coord(image: Molecule, prev_image: Molecule):
        # centering the image
        v = np.mean(np.array([atom.coords for atom in image.pybelobj.atoms]), axis=0)
        image.pybelobj.OBMol.Translate(ob.vector3(-v[0], -v[1], -v[2]))
        image.position = np.array([atom.coords for atom in image.pybelobj.atoms]).flatten()

        # rotation matrix
        c = np.zeros((3, 3))
        for x, y in zip(image.pybelobj.atoms, prev_image.pybelobj.atoms):
            c = c + (np.array(x.coords).reshape(3, 1) @ np.array(y.coords).reshape(1, 3))
        f = np.array(
            [(c[0][0] + c[1][1] + c[2][2], c[1][2] - c[2][1], c[2][0] - c[0][2], c[0][1] - c[1][0]),
             (c[1][2] - c[2][1], c[0][0] - c[1][1] - c[2][2], c[0][1] + c[1][0], c[2][0] + c[2][0]),
             (c[2][0] - c[0][2], c[0][1] + c[1][0], -c[0][0] + c[1][1] - c[2][2], c[1][2] + c[2][1]),
             (c[0][1] - c[1][0], c[2][0] + c[2][0], c[1][2] + c[2][1], -c[0][0] - c[1][1] + c[2][2])])
        eig, v = np.linalg.eig(f)
        q = v[:, eig.argmax()].flatten()
        r = np.array(
            [((q[0] ** 2 + q[1] ** 2 - q[2] ** 2 - q[3] ** 2) / 2, q[1] * q[2] - q[0] * q[3], q[1] * q[3] + q[0] * q[2]),
             (q[1] * q[2] + q[0] * q[3], (q[0] ** 2 - q[1] ** 2 + q[2] ** 2 - q[3] ** 2) / 2, q[2] * q[3] - q[0] * q[1]),
             (q[1] * q[3] - q[0] * q[2], q[2] * q[3] + q[0] * q[1], (q[0] ** 2 - q[1] ** 2 - q[2] ** 2 + q[3] ** 2) / 2)])

        r = 2 * r

        # update coordinate of pybelobj
        for atm in image.pybelobj.atoms:
            rotated_coord = r @ np.array(atm.coords)
            atm.OBAtom.SetVector(rotated_coord[0], rotated_coord[1], rotated_coord[2])
        # update the position
        image.position = np.array([atom.coords for atom in image.pybelobj.atoms]).flatten()

        # store rotation matrix
        r_all = np.array([r for i in range(len(image.pybelobj.atoms))])
        return r_all

    @classmethod
    def rotate_images(cls, images: List[Molecule]):
        matrix = []
        for i in range(1, len(images) - 1):
            rot = cls.rotate_coord(images[i], images[i - 1])
            matrix.append(rot)
        cls.rotate_coord(images[-1], images[-2])
        return matrix


class NEB:
    def __init__(self, spring, method, ci=False):
        self.ksp = spring
        # method should be 'ori', 'reg', 'ci-eb'
        self.method = method
        self.ci = ci

    def nebforce(self, images: List[Molecule]):
        nebforce = []
        if self.method == 'ci-eb':
            leq = np.linalg.norm(images[-1].position - images[0].position, ord=2) / (len(images) - 1)
        else:
            leq = 1

        imax = []
        if self.method == 'ci-eb' or self.ci:
            for k in range(1, len(images) - 1):
                if images[k - 1].energy < images[k].energy and images[k].energy > images[k + 1].energy:
                    imax.append(k)

        for k in range(1, len(images) - 1):
            r_plus = images[k + 1].position - images[k].position
            r_minus = images[k].position - images[k - 1].position
            lr_plus = np.linalg.norm(r_plus, ord=2)
            lr_minus = np.linalg.norm(r_minus, ord=2)

            # define tangent vector
            if self.method == 'reg':
                if images[k - 1].energy < images[k].energy < images[k + 1].energy:
                    tangent = r_plus / lr_plus
                elif images[k - 1].energy > images[k].energy > images[k + 1].energy:
                    tangent = r_minus / lr_minus
                else:
                    v1 = max(abs(images[k + 1].energy - images[k].energy), abs(images[k].energy - images[k - 1].energy))
                    v2 = min(abs(images[k + 1].energy - images[k].energy), abs(images[k].energy - images[k - 1].energy))
                    if images[k - 1].energy < images[k + 1].energy:
                        tangent = v1 * r_plus + v2 * r_minus
                        tangent = tangent / np.linalg.norm(tangent, ord=2)
                    else:
                        tangent = v2 * r_plus + v1 * r_minus
                        tangent = tangent / np.linalg.norm(tangent, ord=2)
            else:
                tangent = r_plus / lr_plus + r_minus / lr_minus
                tangent = tangent / np.linalg.norm(tangent, ord=2)

            # define spring force
            if self.ci and k in imax:
                force_image = images[k].force - 2 * np.dot(images[k].force, tangent) * tangent
            else:
                if self.method == 'reg':
                    force_spring = self.ksp * (lr_plus - lr_minus) * tangent
                elif self.method == 'ci-eb' and (k + 1 in imax or k - 1 in imax):
                    v1 = max(abs(images[k + 1].energy - images[k].energy), abs(images[k].energy - images[k - 1].energy))
                    v2 = min(abs(images[k + 1].energy - images[k].energy), abs(images[k].energy - images[k - 1].energy))
                    force_spring = self.ksp * (
                            (lr_plus - leq) / lr_plus * r_plus - (lr_minus - leq) / lr_minus * r_minus) * v1 / v2
                elif self.method == 'ci-eb':
                    force_spring = self.ksp * (
                            (lr_plus - leq) / lr_plus * r_plus - (lr_minus - leq) / lr_minus * r_minus)
                else:
                    force_spring = self.ksp * np.dot(r_plus - r_minus, tangent) * tangent
                force_pp = images[k].force - np.dot(images[k].force, tangent) * tangent
                force_image = force_spring + force_pp

            nebforce.append(force_image)

        return np.array(nebforce).flatten()


class InitStructure:
    @staticmethod
    def generate_inter_lp(init_str: pybel.Molecule, end_str: pybel.Molecule, inter: int):
        int_list = [pybel.Molecule(ob.OBMol()) for i in range(inter)]
        for i, j in zip(init_str.atoms, end_str.atoms):
            # define the vector for propagation
            i_coord = np.array(i.coords)
            j_coord = np.array(j.coords)
            coord_vector = (j_coord - i_coord) / (inter + 1)
            for k in range(inter):
                # define the atomic position for k-th intermediate and write to OBmol object
                atom_position = i_coord + coord_vector * (k + 1)
                newatom = ob.OBAtom()
                newatom.SetAtomicNum(i.OBAtom.GetAtomicNum())
                newatom.SetVector(atom_position[0], atom_position[1], atom_position[2])
                int_list[k].OBMol.AddAtom(newatom)
        return int_list

    @staticmethod
    def update_idpp_force(mol: Molecule, dist_matrix):
        diff = []
        for i, atom_i in enumerate(mol.pybelobj.atoms):
            force_x = 0.0
            force_y = 0.0
            force_z = 0.0
            for j, atom_j in enumerate(mol.pybelobj.atoms):
                if i != j:
                    d = atom_i.OBAtom.GetDistance(atom_j.OBAtom)
                    df = - 2 / d ** 3 + 6 * dist_matrix[i][j] / d ** 4 - 4 * (dist_matrix[i][j] ** 2) / d ** 5
                    force_x = force_x - df * (atom_i.coords[0] - atom_j.coords[0]) / np.sqrt(d)
                    force_y = force_y - df * (atom_i.coords[1] - atom_j.coords[1]) / np.sqrt(d)
                    force_z = force_z - df * (atom_i.coords[2] - atom_j.coords[2]) / np.sqrt(d)
            diff.extend((force_x, force_y, force_z))
        mol.force = np.array(diff)


class OptimizationMethod(metaclass=ABCMeta):
    @abstractmethod
    def dx_first(self, force: np.ndarray):
        pass

    @abstractmethod
    def dx_opt(self, force: np.ndarray):
        pass


class FIRE(OptimizationMethod):
    # fixed parameters for fire algorithm
    # parameters from PRL 2006, 97, 170201 and dt and dt_max was modified based on a calculation
    dt_max = 1
    n_acc = 5
    f_inc = 1.1
    f_acc = 0.99
    f_dec = 0.5
    a_start = 0.1

    def __init__(self, dim):
        # valuable parameters during optimization
        self.dt = 0.3
        self.a = 0.1
        self.n_prop = 0
        self.v = np.zeros(dim)

    def dx_first(self, force, maxmove=0.2):
        self.v = self.dt * force
        dx = self.dt * self.v
        xmax = np.abs(dx).max()
        if xmax > maxmove:
            step = maxmove / xmax
            dx = step * dx
            self.v = step * self.v
        return dx

    def dx_opt(self, force, maxmove=0.2):
        vp = (1.0 - self.a) * self.v + self.a * np.linalg.norm(self.v) / np.linalg.norm(force) * force
        if np.dot(self.v, force) > 0.0:
            if self.n_prop > self.n_acc:
                self.dt = min(self.dt * self.f_inc, self.dt_max)
                self.a = self.a * self.f_acc
            self.n_prop += 1
        else:
            vp = 0 * vp
            self.a = self.a_start
            self.dt = self.dt * self.f_dec
            self.n_prop = 0
        self.v = vp + self.dt * force
        dx = self.dt * self.v
        xmax = np.abs(dx).max()
        if xmax > maxmove:
            step = maxmove / xmax
            dx = step * dx
            self.v = step * self.v
        return dx


class BFGS(OptimizationMethod):
    def __init__(self, dim):
        self.hesse = np.identity(dim)
        self.stored_forces = np.zeros(dim)
        self.stored_g = np.zeros(dim)

    def dx_first(self, force, maxmove=0.2):
        dx = force
        xmax = np.abs(dx).max()
        if xmax > maxmove:
            step = maxmove / xmax
            dx = step * dx
        self.stored_forces = force
        self.stored_g = dx
        return dx

    def dx_opt(self, force, maxmove=0.2):
        yk = - (force - self.stored_forces)
        sy = np.dot(yk, self.stored_g)
        hesse_new = self.hesse + (1 + yk @ self.hesse @ yk / sy) / sy * np.outer(self.stored_g, self.stored_g)\
                    - (np.outer(self.hesse @ yk, self.stored_g) + np.outer(self.stored_g, yk @ self.hesse)) / sy
        self.hesse = hesse_new
        dx = hesse_new @ force
        xmax = np.abs(dx).max()
        if xmax > maxmove:
            step = maxmove / xmax
            dx = step * dx
        self.stored_forces = force
        self.stored_g = dx
        return dx


class PathOptimization:
    def __init__(self, init_path: List[Molecule], neb: NEB, calcio: GaussianIO, ci=False):
        self.images = init_path
        self.neb = neb
        self.step = 0
        self.rms = 0
        self.calcio = calcio
        self.ci = ci

    def move_structure(self, s_move: np.ndarray, step=1.0):
        for j, image in enumerate(self.images[1:-1]):
            fread_ptr = j * 3 * len(image.pybelobj.atoms)
            for i, atm in enumerate(image.pybelobj.atoms):
                x = atm.coords[0] + step * s_move[fread_ptr + 3 * i]
                y = atm.coords[1] + step * s_move[fread_ptr + 3 * i + 1]
                z = atm.coords[2] + step * s_move[fread_ptr + 3 * i + 2]
                atm.OBAtom.SetVector(x, y, z)
            # update coord
            image.position = np.array([atom.coords for atom in image.pybelobj.atoms]).flatten()

    def write_energy(self, file):
        # change directory and write energy to result.csv
        os.chdir('..')
        with open(file, mode='a') as fr:
            fr.write(
                '{0},{1},{2}\n'.format(str(self.step), ','.join([str(i.energy) for i in self.images]), str(self.rms)))

        # change to working directory
        os.chdir('temp')

    def write_molecules(self):
        # change directory and create xyz file
        os.chdir('..')
        outstream = pybel.Outputfile("xyz", "path" + str(self.step) + ".xyz")
        for mol in self.images:
            mol.pybelobj.OBMol.SetTitle(mol.name + " Energy=" + str(mol.energy))
            outstream.write(mol.pybelobj)
        outstream.close()

        # change to working directory
        os.chdir('temp')

    def optimize(self, maxiter, threshold, method: OptimizationMethod):
        # change directory and initialize result file
        os.chdir('..')
        result_file = 'result.csv'
        with open(result_file, mode='w') as fr:
            fr.write('step,start,{0},end,rms_force\n'.format(
                ','.join(['int_' + str(i + 1) for i in range(len(self.images) - 2)])))

        # change to working directory
        os.chdir('temp')

        while self.step != maxiter:
            # force calculation
            for mol in self.images[1:-1]:
                self.calcio.update_force(mol, 'Force') if self.step == 0 \
                    else self.calcio.update_force(mol, 'Force guess=read')
            nebforce = self.neb.nebforce(self.images)
            self.rms = np.sqrt(np.mean(nebforce ** 2))

            # write path N
            self.write_energy(result_file)
            self.write_molecules()

            # check conv.
            if self.rms < threshold:
                return True

            # if ci-neb, change from regular neb after 20 steps
            if self.ci and self.step == 20:
                self.neb.ci = True

            # move structure
            dx = method.dx_first(nebforce) if self.step == 0 else method.dx_opt(nebforce)
            self.move_structure(dx)
            self.step += 1

        return False

    def optimize_idpp(self):
        # value
        maxiter = 200
        threshold = 0.02
        dim = (len(self.images) - 2) * 3 * len(self.images[0].pybelobj.atoms)
        method = FIRE(dim)

        # calculate distance lists
        dm_init = np.array([[i.OBAtom.GetDistance(j.OBAtom) for j in self.images[0].pybelobj.atoms]
                            for i in self.images[0].pybelobj.atoms])
        dm_end = np.array([[i.OBAtom.GetDistance(j.OBAtom) for j in self.images[-1].pybelobj.atoms]
                           for i in self.images[-1].pybelobj.atoms])
        dm_diff = (dm_end - dm_init) / (len(self.images) - 1)
        dm_int = [dm_init + (i + 1) * dm_diff for i in range(len(self.images) - 2)]

        while self.step != maxiter:
            for mol, dm in zip(self.images[1:-1], dm_int):
                InitStructure.update_idpp_force(mol, dm)
            nebforce = self.neb.nebforce(self.images)
            self.rms = np.sqrt(np.mean(nebforce ** 2))
            print(str(self.rms))

            # check conv.
            if self.rms < threshold:
                break
            else:
                dx = method.dx_first(nebforce, 0.1) if self.step == 0 else method.dx_opt(nebforce, 0.1)
                self.move_structure(dx)
                self.step += 1

        return False


def main(file):
    # extract input information from gjf file
    with open(file) as f1:
        strs = f1.read()

    reg = re.compile(
        r'%chk=(.*).chk\n(.*)(#.+?)\n *\n.*?(|!neb=(.+?))\n *\n(.+?)\n(.+?\n) *\n.*?\n *\n.+?\n(.+?\n) *\n(.*)',
        re.MULTILINE | re.DOTALL | re.IGNORECASE)
    match = reg.search(strs)
    if match:
        neb_opt = match.group(5)
        cm = match.group(6)
        init_coord = match.group(7)
        end_coord = match.group(8)
        calc = CalcMethods(match.group(1), match.group(2), match.group(3), match.group(9))
    else:
        print("cannot read input file")
        return -1

    # Gaussian Parser
    calcio = GaussianIO(calc)

    # default values
    inter = 10
    max_iter = 20
    spring = 0.02  # unit: Hartree / Bohr^2
    conv_threshold = 0.001
    neb = 'reg'
    init = 'linear'
    ci = False
    opt = 'bfgs'

    # read NEB option
    options = re.findall(
        r'inter=(\d+)|maxiter=(\d+)|spring=(\d+)|conv=(\d+)|neb=(reg|org|ci-eb)|init=(idpp|linear)|(ci)|opt=(BFGS|FIRE)',
        neb_opt, re.IGNORECASE)
    if options:
        for option in options:
            if option[0] != '':
                inter = int(option[0])
            if option[1] != '':
                max_iter = int(option[1])
            if option[2] != '':
                spring = float(option[2])
            if option[3] != '':
                conv_threshold = float(option[4])
            if option[4] != '':
                neb = option[4].lower()
            if option[5] != '':
                init = option[5].lower()
            if option[6] != '':
                ci = True
            if option[7] != '':
                opt = option[7].lower()

    # NEB method
    nebmethod = NEB(spring, neb)

    # make xyz-file input stream
    numatom = init_coord.count('\n')
    mol_init = Molecule(pybel.readstring("xyz", str(numatom) + "\n" + "init\n" + init_coord), 'start', cm)
    mol_end = Molecule(pybel.readstring("xyz", str(numatom) + "\n" + "end\n" + end_coord), 'end', cm)

    # centering of start coordinate
    v = np.mean(np.array([atom.coords for atom in mol_init.pybelobj.atoms]), axis=0)
    mol_init.pybelobj.OBMol.Translate(ob.vector3(-v[0], -v[1], -v[2]))
    mol_init.position = np.array([atom.coords for atom in mol_init.pybelobj.atoms]).flatten()
    # centering and rotation of end coordinate
    Rotate.rotate_coord(mol_end, mol_init)

    # Generate initial intermediate structures and store them into List
    mol_beads = [mol_init]
    for k, mol in enumerate(InitStructure.generate_inter_lp(mol_init.pybelobj, mol_end.pybelobj, inter)):
        mol_beads.append(Molecule(mol, 'int_' + str(k + 1), cm))
    mol_beads.append(mol_end)

    # initial structure optimization by idpp
    if init == 'idpp':
        PathOptimization(mol_beads, NEB(0.02, 'org'), calcio).optimize_idpp()

    # create working directory and go into it
    os.mkdir('temp')
    os.chdir('temp')

    # run Gaussian calculation for start and end point
    calcio.update_energy(mol_beads[0])
    calcio.update_energy(mol_beads[-1])

    # optimization step
    dim = (len(mol_beads) - 2) * 3 * len(mol_beads[0].pybelobj.atoms)
    if opt == 'fire':
        opt_flag = PathOptimization(mol_beads, nebmethod, calcio, ci).optimize(max_iter, conv_threshold, FIRE(dim))
    else:
        opt_flag = PathOptimization(mol_beads, nebmethod, calcio, ci).optimize(max_iter, conv_threshold, BFGS(dim))

    if opt_flag:
        print("Result: Converged")
    else:
        print("Result: Not converged")


if __name__ == '__main__':
    import sys

    args = sys.argv
    # check the input
    if len(args) != 2:
        print("usage: nev-script.py <input.gjf>")
    main(args[1])
