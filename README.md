# OSMES
Pipeline in python to perform a One Substrate-Many Enzymes Screening (OSMES) of substrates bound as external aldimine towards PLPomes
<br><br>
![OSMES_H](./OSMES_H.png)
<br>
## Installation 
### Requirements
Any version of [mamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)

### Cloning git repository
```bash
git clone https://github.com/lab83bio/OSMES.git
cd OSMES
```
### Installation of `ADFRsuite_x86_64Linux_1.096`
```bash
wget https://ccsb.scripps.edu/adfr/download/1038/ -O 'adfr.tar.gz'
tar zxvf adfr.tar.gz 
cd ADFRsuite_x86_64Linux_1.0/
chmod +x install.sh
yes|bash install.sh -d ../ADFRsuite-1.0 -c 0 &>../ADFR_install.log
cd ..
```
### Create and activate mamba envinroment
```bash
micromamba create -f OSMES_explicit_env.txt -n OSMES -y
micromamba activate OSMES
```
## Usage
### Make input files
OSMES pipleine required the following input files:
- `substrate.pdbqt` for ADFR
- `substrate.txt` for calculation of the catalytic favourable conformations (CFC)
- Enzyme set in pdb format
- `coords.tsv` with the specificied coordinates of the gridbox center for each active site.
  
To produce the indicated input files, please use the [*OSMES.ipynb*](https://colab.research.google.com/drive/1lF4ezjLnJ16w6RrC5R_5ZV0P5g9omAtd#scrollTo=AfUiKQWES7V8) Colab Notebook and follow the instruction of the cells.
### Usage of OSMES
make executable the file
```bash
chmod +x OSMES_submit.py
```
and then run the pipeline with the configuration file for [test](https://github.com/lab83bio/OSMES/tree/main/test)
```bash
OSMES_submit.py OSMES.config
```
where in OSMES config are defined all the path and the parameters required by the analysis, an example below:
```python
[PATHS]
# Specific input/output files
adfr_path = ADFRsuite-1.0/bin
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
