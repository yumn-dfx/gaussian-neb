import re
from openbabel import openbabel as ob
from openbabel import pybel
import nebscript as ns
import numpy as np


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
        calc = ns.CalcMethods(match.group(1), match.group(2), match.group(3), match.group(9))
    else:
        print("cannot read input file")
        return -1

    # Gaussian Parser
    calcio = ns.GaussianIO(calc)

    # default values
    inter=10
    init = 'linear'

    # read NEB option
    options = re.findall(r'inter=(\d+)|init=(ldpp|linear)|', neb_opt, re.IGNORECASE)
    if options:
        for option in options:
            if option[0] != '':
                inter = int(option[0])
            if option[1] != '':
                init = option[1].lower()

    # make xyz-file input stream
    numatom = init_coord.count('\n')
    mol_init = ns.Molecule(pybel.readstring("xyz", str(numatom) + "\n" + "init\n" + init_coord), 'start', cm)
    mol_end = ns.Molecule(pybel.readstring("xyz", str(numatom) + "\n" + "end\n" + end_coord), 'end', cm)

    # centering of start coordinate
    v = np.mean(np.array([atom.coords for atom in mol_init.pybelobj.atoms]), axis=0)
    mol_init.pybelobj.OBMol.Translate(ob.vector3(-v[0], -v[1], -v[2]))
    mol_init.position = np.array([atom.coords for atom in mol_init.pybelobj.atoms]).flatten()
    # centering and rotation of end coordinate
    ns.Rotate.rotate_coord(mol_end, mol_init)

    # Generate initial intermediate structures and store them into List
    mol_beads = [mol_init]
    for k, mol in enumerate(ns.InitStructure.generate_inter_lp(mol_init.pybelobj, mol_end.pybelobj, inter)):
        mol_beads.append(ns.Molecule(mol, 'int_' + str(k + 1), cm))
    mol_beads.append(mol_end)

    # ldpp initial structure generation
    if init == 'ldpp':
        ns.PathOptimization(mol_beads, ns.NEB(0.02, 'org'), calcio).optimize_ldpp()

    # write the gjf file
    calcio.write_com('trial', mol_beads[1].pybelobj, mol_beads[1].cm)

    # write initial trj into xyz file
    outstream = pybel.Outputfile("xyz", "initpath.xyz", overwrite=True)
    for mol in mol_beads:
        mol.pybelobj.OBMol.SetTitle("init_str")
        outstream.write(mol.pybelobj)
    outstream.close()

    return 0


if __name__ == '__main__':
    import sys

    args = sys.argv
    # check the input
    if len(args) != 2:
        print("usage: init_str.py <input.gjf>")
    main(args[1])
