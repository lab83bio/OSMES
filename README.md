# RevDockPLP
Pipeline in python to perform a reverse docking screening with PLPomes and external aldimine as ligands

## Requirements
Any version of micromamba or mamba https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html

## Local utilizaiton 
### Cloning git repositories
```bash
mkdir OSMES
cd OSMES
git clone https://github.com/lab83bio/RevDockPLP.git
```
### Installing `ADFRsuite_x86_64Linux_1.096`
```bash
wget https://ccsb.scripps.edu/adfr/download/1038/ -O 'adfr.tar.gz'
tar zxvf adfr.tar.gz 
cd ADFRsuite_x86_64Linux_1.0/; chmod +x install.sh; yes|./install.sh -d ../ADFRsuite-1.0 -c 0 &>../ADFR_install.log
cp RevDockPLP/ade.py ADFRsuite_x86_64Linux_1.0/ADFRsuite-1.0/CCSBpckgs/ADFR/bin
```
### Create and activate mamba envinroment
```bash
micromamba create -f OSMES_explicit_env.txt -n OSMES -y
micromamba activate OSMES
```
## Creating input files
To create input files such as ligand or enzyme set, please use the Colab Notebook below:
https://colab.research.google.com/drive/1lF4ezjLnJ16w6RrC5R_5ZV0P5g9omAtd#scrollTo=AfUiKQWES7V8

## Usage of OSMES
```bash
OSMES_submit.py OSMES.config
```
where in OSMES config are defined all the path and the parameters required by the analysis, an example below:
```python
[PATHS]
# Specific input/output files
adfr_path = ../ADFRsuite-1.0/bin
Ligand = test/HTL_PLP.pdbqt # substrate file
reactions_file = test/HTL_PLP.txt # file with the atoms for the catalyitic favourable conformations and the gridbox sizes
receptor_dir = test/Mus_musculus/ # path to enzyme pdb files
coord_file = test/Mus_musculus/Mus_musuclus_coord.csv # coordinates file for every active site consdered in the enzymes set with specified center of the gridbox
outdir = test/OSMES_results # directory of results

[GPF]
spacing=0.375 
pocketMode = all

[DPF]
maxEvals = 5000000
nbRuns = 20
maxCores = 0 # Set 0 to use all cores available
clusteringRMSDCutoff = 3
popSize = 300 

[OTHER]
imgFormat = pdf
flex_lysine = True # keep catalytic lysine as flexible receptor
lig_res = HTL # substrate residue name of pdbqt file
```
