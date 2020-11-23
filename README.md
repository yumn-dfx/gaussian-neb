A small script to do NEB calculation by Gaussian (alpha ver.) 
=======================================

What's this?
-------------------------
pythonの勉強を兼ねて書いてみたスクリプト，試作品  
どうせ書くならありそうでないものを書いてみようかなと  
とはいえNWChem，turbomole，ORCAなどのソフトがNEB計算に対応しますし，GaussianでやりたいならReactionPlusProもあるのであまり存在価値がないような


Requirement
------------------------
おそらくpython 3.5-   
numpy  
openbabel  
pybel (openbabelと一緒にインストールされるはず)  
Gaussian09(おそらく16も可)  
  
動作確認は以下の環境で行っています  
(miniconda環境)  
python 3.8  
openbabel 3.1  
numpy 1.19  
Gaussian09 revD01  
  
Gaussian16環境がないのでGaussian16での動作確認はしていませんが，output fileの形式は変わっていないと思うので，おそらく動くと思われます．

Usage
------------------------
Gaussian16を使う場合30行目あたりにある
~~~
class GaussianIO:
    # change the following value to 'g16' if you want to run with Gaussian 16
    gau_cmd = 'g09'
~~~
を書き換えてください．  
また，openbabelのversionが2.4以下の場合，最初の
~~~
from openbabel import openbabel as ob
from openbabel import pybel
~~~
を
~~~
import openbabel as ob
import pybel
~~~
に書き換えてください．  

gjf(あるいはcom)ファイルのあるディレクトリに移動し  
~~~
hogehoge$ python <path to the script>/nebscript.py input.gjf
~~~
カレントディレクトリにtempフォルダを作成し，その中に中間構造の計算のためのgjfファイルおよび出力のoutファイルを生成します．  

補助スクリプトとして，最初の中間構造を出力するinit_str.pyがあります．このスクリプトを実行すると読み込んだ構造を
もとに初期の中間構造をinitpath.xyzに，gjfファイルの出力例としてtrial.gjfを書き出します．  
初期の中間構造にあまりに無理のあるとHF/DFT計算ができませんので，補助スクリプトを使って無理のない中間体の構造が出力
されているか確認されるとよいかと思います．

### gjfファイルに関して

基本的にはGaussianでのQST2計算のInputファイルと同じですが，かなり適当なgjfファイルの読み込みをしている関係上すこし制限があります．  
具体的に，Gaussianでは%chkコマンドが１行目にある必要はありませんが，読み込みの関係上１行目に指定してください．  
現地点では，始点と終点の２つの構造の読み込みしか対応していません．また，conectivityの読み込みも対応していません．  
>gjfファイルの例  
~~~
%chk=(..../)<<filename>>.chk   <-必須
%mem=1000MB                    <-使用メモリー量などのLink 0コマンド:オプション
......  
# b3lyp/6-31G ....             <-計算手法　注：opt, freq等のキーワードを入れないで下さい
(空行)
title !neb=inter=10　           <-titleのあとにNEB計算のオプションを指定：オプション
(空行)
0 1                            <- 始点となる構造の構造情報
C 1.00000 1.00000 1.000000   
.
.
N 5.00000 5.00000 5.000000
(空行)
title
(空行)
0 1                            <- 終点となる構造の構造情報　ここで指定された電荷と多重度は無視しています
C 2.00000 2.00000 2.000000   
.
.
N 4.00000 4.00000 4.000000
(空行)
P 0                            <- 基底などについて追加の指定：必要に応じて
6-31+G* 
~~~

### NEB計算のオプション
NEB計算のオプションを指定する際は以下のように指定して下さい．かっこはあってもなくても．  
指定がなければデフォルト値を使用します．!nebの記述がなければ，すべてデフォルト値で計算を始めます．  
~~~
!neb=(inter=8,maxiter=30,spring=0.02,step=1,conv=0.0001,init=idpp,ci)  
~~~
#### 各オプションについて  
- inter=8  始点と終点を結ぶ間の構造の数　デフォルト:10  
- maxiter=20  NEB最適化をする回数の最大値　デフォルト:20  
- spring=0.02 構造間を結ぶバネのばね定数(単位 Hartree/Bohr**2)　デフォルト:0.02  
- conv=0.0001 収束判定条件:各構造にかかるForceのRMSのしきい値　デフォルト:0.001  
- init=idpp 初期構造の生成方法  
linear: 始点と終点の座標から直線的に初期構造を生成します．デフォルト  
idpp: idpp法(J. Chem. Phys. 2014, 140, 214106.)による初期構造最適化をします(試験的)．  
linearで無茶な構造が出力されたときに使うことをおすすめします．
- neb=reg 中間構造にかかる力の計算方法  
org: オリジナルの方法で力を計算します  
reg: 改良された方法(J. Chem. Phys. 2000, 113, 9978.)で力を計算します．デフォルト  
ci-eb: CI-EB法(J. Chem. Phys. 2016, 145, 094107.)で力の計算をします（下記のCIオプションは不要です）  
- opt=fire　最適化の方法  
fire: FIRE法で最適化します．
bfgs: BFGS法で最適化します．デフォルト
- ci　このオプションが指定された場合，Climbing Image NEBを行います  


Output
------------------------
tempフォルダ内に生成したGaussianのinputファイルおよびGaussianからの出力を書き出します．  
各最適化ステップにおいて，中間体の構造をpath\<step number>.xyzファイルに出力します．  
また，result.csvに各最適化ステップでの始点，中間構造，終点のエネルギー値および力のRMS値を出力します．  

Comments
------------------------
勉強がてら書いてみただけの簡単なスクリプトで著作権に関してどうこう言うつもりもないですが，MIT Licenseにしておきます．  
ただし，openbabelがGPLなので，MITにしていいのか自信がありません．(詳しいかたがいましたらご指摘いただければ)  
最近は量子化学計算をすることもあまりなくなったので，Issueを上げていただいてもあまり対応できないかもしれません．  