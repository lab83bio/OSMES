#!/usr/bin/env python3
import configparser, sys, os

#####    config parser    ####
config = configparser.ConfigParser(inline_comment_prefixes="#")
config.optionxform = str
try:
    config.read(sys.argv[1])
    locals().update(dict(config.items('PATHS')))
    locals().update(dict(config.items('GPF')))
    locals().update(dict(config.items('DPF')))
    locals().update(dict(config.items('OTHER')))    
except:
    print(
'''
### Configuration file for the reverse docking procedure ###
###
### change file names according with your paths ####
###

[PATHS]
# Specific 'input/output files'
adfr_path = f'/hpc/home/marco.malatesta/revdocking/ADFRsuite_x86_64Linux_1.0/ADFRsuite-1.0/bin'
Ligand = /hpc/home/marco.malatesta/Rev_Docking/Docking_data/ORN_PLP.pdbqt            # ligand file
reactions_file = 'HTL_PLP.txt'
receptor_dir = /hpc/archive/G_BIOSCIENZE/marco.malatesta/PLP_enzymes/Mus_musculus/AF_db/ # path to receptor pdb files
coord_file = /hpc/archive/G_BIOSCIENZE/marco.malatesta/PLP_enzymes/Mus_musculus/AF_db/Mus_musculus_coord.csv  outdir = OMSES_results # directory of results
outdir = test/OSMES_results # directory of results

[GPF]
npts=20 20 20
spacing=0.2
pocketMode = "all"

[DPF]
maxEvals = 50000000
nbRuns = 200
maxCores = 0 # Set 0 to use all cores available
clusteringRMSDCutoff = 3
popSize = 300 

[OTHER]
imgFormat = 'pdf'
flex_lysine = True 
lig_res = 'HTL'
'''
    )
    sys.exit(1)

### ABSOLUTE PATHS
adfr_path = os.path.abspath(adfr_path)
Ligand = os.path.abspath(Ligand)
reactions_file = os.path.abspath(reactions_file)
receptor_dir = os.path.abspath(receptor_dir)
coord_file = os.path.abspath(coord_file)
rec_files_dir = os.path.abspath(receptor_dir)
outdir = os.path.abspath(outdir)

from OSMES import *
warnings.filterwarnings('ignore')
cwd = os.getcwd()
os.chdir(cwd)
rec_files_dir = ''
## Setting varaiables
lig_name = os.path.splitext(os.path.basename(Ligand))
dihedral = eval(open(reactions_file).read())

adfr_help = os.popen(f'{adfr_path}/adfr -h').readlines()
agfr_help = os.popen(f'{adfr_path}/agfr -h').readlines()
adfr_flags = [x.split('--')[1].split(' ')[0] for x in adfr_help if x.strip().startswith('-')]
agfr_flags = [x.split('--')[1].split(' ')[0] for x in agfr_help if x.strip().startswith('-')]
adfr_options = ' '.join([f"--{i} {a}" for i, a in locals().items() if i in adfr_flags if i != 'config'])
agfr_options = ' '.join([f"--{i} {a}"  for i, a in locals().items() if i in agfr_flags if i != 'config'])
print(agfr_options, adfr_options)
size = max([14, dihedral[4]]) #@param {type:"slider", min:10, max:60, step:5} 
npts = ' '.join([str(size)]*3) # npts = "25 25 25"

i = 1
today = date.today()
outdir_date = f"{outdir}_{today.strftime('%d-%m-%Y')}"

while os.path.exists(outdir_date+'_%s' % i):
    i += 1
out_dir = f"{outdir_date}_{i}/"
os.system(f'mkdir -p {out_dir}')
if os.path.exists(outdir_date+'_%s' % i):
    print('writing output in '+ out_dir)

####    read coordinates    ####
print(f'Reading coordinates from {coord_file} ...')
df_coord = pd.read_table(coord_file, comment='#').dropna()
# df_coord = df_coord[df_coord['pdb'].str.contains('F222')]

####   copy the Ligand   ####
os.system(f'cp {Ligand} {out_dir}')
lig_name = os.path.splitext(os.path.basename(Ligand))[0]


####    reverse docking    ####
for i, r in df_coord.iterrows():
    pdb = r['pdb']
    xyz = r[['x', 'y', 'z']]
    pdb_name = os.path.splitext(pdb)[0]
    lys = r['res']
    i = r['chain'] 
    try:
        trg_file = f'{pdb_name}_{i}.trg'
        print(f"## Docking {lig_name} in **{pdb_name}** using LYS {int(lys)} in chain {i}")
        if rec_files_dir:
            # os.chdir(cwd)
            os.system(f'cp {rec_files_dir}/{pdb_name}_{i}.trg {out_dir}')
            os.system(f'cp {rec_files_dir}/{pdb_name}_{i}.pdbqt {out_dir}')
            # os.chdir(out_dir)
        else:

            ####    prepare receptor    ####
            os.chdir(cwd)
            print(f'### Prepare {pdb_name}_{i}.pdbqt ...')
            print(re.sub('HETATM', '#HETATM', open(f'{receptor_dir}/{pdb}', 'r').read()), 
                file=open(out_dir+pdb, 'w'))                   # comment HETATM in pdb file
            os.chdir(out_dir)
            subprocess.run([f'{adfr_path}/prepare_receptor',
                            '-A', 'checkhydrogens',
                            '-e',
                            '-r', f'{receptor_dir}/{pdb}',
                            '-o', f'{out_dir}/{pdb_name}_{i}.pdbqt'
                            ],
                        stdout=open(f"{pdb_name}_{i}_prepare_receptor.log", 'w'),
                        stderr=open(f"{pdb_name}_{i}_prepare_receptor.err", 'w')
                        )
            if os.stat(f"{pdb_name}_{i}_prepare_receptor.err").st_size != 0:
                continue
            os.chdir(cwd)

            ####    AGFR    ####
            print(f'### Run AGFR for grid computation ...')
            os.chdir(out_dir)
            if flex_lysine:
                subprocess.run([f'{adfr_path}/agfr',
                                '-r', f'{pdb_name}_{i}.pdbqt',
                                '-b', f'user {" ".join(map(lambda x: str(round(x,3)),xyz))} {npts} {agfr_options}',
                                '-f', f'{i}:LYS{str(int(lys))}'
                                ],
                            stdout=open(f"{pdb_name}_{i}_agfr.log", 'w'),
                            stderr=open(f"{pdb_name}_{i}_agfr.err", 'w')
                            )
            else:
                subprocess.run([f'{adfr_path}/agfr',
                                '-r', f'{pdb_name}_{i}.pdbqt',
                                '-b', f'user {" ".join(map(lambda x: str(round(x,3)),xyz))} {npts} {agfr_options}',
                                ],
                            stdout=open(f"{pdb_name}_{i}_agfr.log", 'w'),
                            stderr=open(f"{pdb_name}_{i}_agfr.err", 'w')
                            )
            if os.stat(f"{pdb_name}_{i}_agfr.err").st_size != 0:
                continue
            os.chdir(cwd)

        ####    ADFR    ####
        print(f'### Run ADFR ...')
        os.chdir(out_dir)
        subprocess.run([f'{adfr_path}/adfr',
                        '-t', f'{pdb_name}_{i}.trg',
                        '-l', f'{os.path.basename(Ligand)}',
                        '-J', f'{pdb_name}_{i}',
                        adfr_options,
                        '--fullOutput',
                        ],
                    stdout=open(f"{pdb_name}_{i}_adfr.log", 'w'),
                    stderr=open(f"{pdb_name}_{i}_adfr.err", 'w')
                    )
        if os.stat(f"{pdb_name}_{i}_adfr.err").st_size != 0:
            continue
        os.chdir(cwd)

        #### PLOT ####
        os.chdir(out_dir)
        conformations = glob.glob(f'*{pdb_name}_{i}*//*[0-9].pdbqt')
        df = pd.concat([calc_run(x, dihedral) for x in conformations])

        # compute ade energies
        ade = pd.DataFrame([calc_ade(c, trg_file) for c in conformations], columns =['run_file', 'ade']) 

        # ade = pd.DataFrame({'energy': [], 'run_file':[]})
        df = df.merge(ade, on='run_file', how='outer')

        histo = DLGdf(glob.glob(f'*{pdb_name}_{i}*summary.dlg')[0], df, lig_res)
        histo.to_csv(f"{pdb_name}_{i}.tsv", sep ='\t', index=False)
        ########

        reactions = [f'{chr(948)} (A)',f'{chr(948)} (B)', f'{chr(948)} (D)']
        reactions = ['X1','X2', 'X3']
        CFCs = ['CFC_1', 'CFC_2', 'CFC_3']
        CFCsBool = ['CFC_1_bool', 'CFC_2_bool', 'CFC_3_bool']

        histo[CFCsBool] = histo[['d', *reactions]].apply(lambda x: x[reactions]==max(x[reactions]) 
                if all([
                     # sorted(x[reactions])[-1] - sorted(x[reactions])[-2] >=0.1, 
                         x['d']<=5
                        ]) 
                else x[reactions]==None, axis=1) 
        histo_mean = pd.concat(list(
                itertools.chain(*[[histo.groupby('rank')[[*reactions, 'd']].mean(),
                                    histo.groupby('rank')['d'].count()],
                                    [pd.DataFrame(histo.groupby('rank')[c].apply(lambda x: len(x[x]))) for c in histo[CFCsBool]]
                                    ])), axis=1)
        histo_mean.columns = ['X1_mean', 'X2_mean', 'X3_mean', 'd mean', 'runs', *CFCs]
        histo = histo_mean.reset_index().merge(histo[['file', 'rank', 'LCC', 'affinity', 
                                                # 'ade_min', 'ade_mean'
                                               ]], on='rank')
        fig, ax = plotCatalyitic(histo)
        fig.savefig(f'catalytic_plot_{pdb_name}_{i}.{imgFormat}', dpi=600)
        os.chdir(cwd)
        ###############

    #     histo = histo[[*CFCs, 'rank', 'runs', 'affinity']]
    #     histo['cat_runs'] = histo[CFCs].apply(list, axis=1)
    #     histo['reaction'] = [['decarboxylase', 'betayase', 'aminotransferase']]*len(histo)
    #     histo = histo.explode(['reaction', 'cat_runs'])
    #     sns.barplot(data=histo, x='rank', y='runs', color='yellow', linewidth=1, edgecolor=".2", ax=ax)
    #     sns.barplot(data=histo, x='rank', y='cat_runs', color='red', hue='reaction', linewidth=1, edgecolor=".2", ax=ax)
        # os.chdir(cwd)
    except:
        print(f'{pdb}_{i} ERROR')
        os.chdir(cwd)

            
            
            
            
            
            
            
            
