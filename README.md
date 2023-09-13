# OSMES
Pipeline in python to perform One Substrate-Many Enzymes Screening (OSMES) of a substrate molecule bound to a PLP cofactor (external aldimine) against a library of PLP-dependent enzyme structures
<br><br>
![OSMES_H](./OSMES_H.png)
<br>
## Installation 
### Requirements
#### OS Requirements
The package development version is tested on Linux operating systems.
#### Software Requirements
Any version of [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) >1.0 <br>
Packages requirement are listed in the [env file](https://github.com/lab83bio/OSMES/blob/main/OSMES_explicit_env.txt)
#### Hardware Requirements
>10 cores
>10 Gb RAM
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
Installation time: ~5 min
### Create and activate mamba environment
```bash
micromamba create -f OSMES_explicit_env.txt -n OSMES -y
micromamba activate OSMES
```
## Usage
### Input files
OSMES pipeline required the following input files:
- `substrate.pdbqt` for ADFR
- `substrate.txt` for calculation of the catalytic favourable conformations (CFC)
- Enzyme set in pdb format
- `coords.tsv` with the specified coordinates of the gridbox center for each active site.
  
To produce the indicated input files, please use the [*OSMES.ipynb*](https://colab.research.google.com/drive/1lF4ezjLnJ16w6RrC5R_5ZV0P5g9omAtd#scrollTo=AfUiKQWES7V8) Colab Notebook and follow the instructions of the cells.
### Running OSMES
make executable the file
```bash
chmod +x OSMES_submit.py
```
and then run the pipeline with the configuration file using the [test dataset](https://github.com/lab83bio/OSMES/tree/main/test).
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
flex_lysine = True # keep catalytic lysine as flexible receptor.
```
Expected output for the test set is in [test/OSMES_results](https://github.com/lab83bio/OSMES/tree/main/test/OSMES_results) directeory
Expected runtime ~30 min 
