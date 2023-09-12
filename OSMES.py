# /usr/bin/env python3
import pandas as pd, gzip, tarfile, urllib, os, sys, glob, requests, json, random, numpy as np, re, time, itertools
from datetime import date

from IPython.display import Markdown, display

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers

from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Blast import NCBIXML
from Bio.PDB import PDBParser as PDB, PDBIO
from Bio.PDB.vectors import *
from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from pymol import cmd, stored
from biopandas.pdb import PandasPdb
import py3Dmol

import pubchempy 
from io import StringIO
from multiprocessing import Pool
from subprocess import Popen, PIPE
import subprocess
from math import sqrt

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

import warnings

def printmd(string):
    display(Markdown(string))
    
def mol_with_atom_index(mol):
    atoms = mol.GetNumAtoms()
    for idx in range( atoms ):
        mol.GetAtomWithIdx( idx ).SetProp( 'molAtomMapNumber', str( mol.GetAtomWithIdx( idx ).GetIdx() ) )
    return mol

def get_structure(mol,n_confs):
    mol = Chem.AddHs(mol)
    new_mol = Chem.Mol(mol)

    AllChem.EmbedMultipleConfs(mol,numConfs=n_confs,useExpTorsionAnglePrefs=True,useBasicKnowledge=True, numThreads=0)
    energies = AllChem.MMFFOptimizeMoleculeConfs(mol,maxIters=2000, nonBondedThresh=100.0)

    energies_list = [e[1] for e in energies]
    min_e_index = energies_list.index(min(energies_list))

    new_mol.AddConformer(mol.GetConformer(min_e_index))

    return new_mol

def convert_ac(list, From, To, db = None):
    '''
    IN: list of accession number (NCBI,UniProt; ENSEMBL, orthodb, ecc)
    original DB (see: https://www.uniprot.org/help/api_idmapping)
    new DB (see: https://www.uniprot.org/help/api_idmapping)
    taxid (eg. 9606), DEFAULT=all
    OUT: dataframe with ID mapping of the 2 database selected
    '''
    url = 'https://rest.uniprot.org/idmapping/run'

    params = {
        'from': From,
        'to': To,
        'ids': ','.join(list),
        'organism_id': str(taxid)
        }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    jobid = json.load(urllib.request.urlopen(req)).get('jobId')
    df = pd.read_table(f"https://rest.uniprot.org/idmapping/results/{jobid}?format=tsv&size=500")
    if db:
        df = pd.read_table(f"https://rest.uniprot.org/idmapping/{db}/results/{jobid}?format=tsv&size=500&")

    return df

def pdb2fasta(pdb):
    """convert pdb file to fasta format for each chain"""
    fasta = []
    
    if pdb.startswith('http'):
        if pdb.endswith('.gz'):
            seq = list(SeqIO.parse(gzip.open(urllib.request.urlopen(pdb),'rt'), 'pdb-atom'))
        else:
            seq = list(SeqIO.parse(StringIO(urllib.request.urlopen(pdb).read().decode('UTF-8')), 'pdb-atom'))

    elif pdb.endswith('.gz'):
        seq = list(SeqIO.parse(gzip.open(pdb,'rt'), 'pdb-atom'))
        
    else:
        seq = list(SeqIO.parse(open(pdb,'r'), 'pdb-atom'))
        
    for record in seq:
        
        if record.annotations['start'] < 1:
            fasta.append(SeqRecord((record.seq[abs(record.annotations['start'])+1:]), 
                             description = pdb.strip('.pdb'), 
                             id = record.annotations['chain']))
        else:
            fasta.append(SeqRecord(('X'*int(record.annotations['start']-1)+record.seq[0:]),
                             description = pdb.strip('.pdb'),
                             id = record.annotations['chain']))
    return fasta


def match_fasta_position(query,subject_list,num=None):
    """return dataframe of position matching in subject and query fasta, or a dictionary given a list of number"""
    
    subject = '\n'.join([f'>{x.id}\n{str(x.seq)}' for x in subject_list])
    child = Popen(' '.join(['blastp', 
                         '-query',fr'<(echo -e ">{query.id}\n{query.seq}")',
                         '-subject', fr'<(echo -e "{subject}")',
                         '-evalue', '10e-5',
                         '-outfmt', '5',
                         '-max_hsps', '1']),
                        stdout=PIPE, stderr=PIPE, shell=True, executable='/bin/bash'
           ) 
    
    df = []
    xml = list(SearchIO.parse(child.stdout, "blast-xml"))
    for q in xml:
        for s in q:
            x=s[0]
            hit_gap = np.array([index for index, value in enumerate(x.hit) if value == '-'])
            q_gap = np.array([index for index, value in enumerate(x.query) if value == '-'])
            hit_num=list(np.arange(x.hit_start+1,x.hit_end+1))
            q_num=list(np.arange(x.query_start+1,x.query_end+1))
            for i in hit_gap:
                hit_num.insert(i,np.nan)
            for i in q_gap:
                q_num.insert(i,np.nan)
            df.append(list(zip(len(x.query)*[query.id],len(x.query)*[s.description],x.query,q_num,x.hit,hit_num)))
    df = pd.DataFrame([j for i in df for j in i],columns=['query','sequence','query_res','query_num','hit_res','hit_num'])
    df.set_index(['query'],inplace=True)
    if num==None:
        return df
    else:
        return df[df.query_num.isin(list(map(int,num)))].to_dict('records')
    

def build_homoligo(subject, template, ac='structure'):
    cmd.reinitialize()
    try:
        cmd.load(template, 'template')
    except:
        cmd.read_pdbstr(urllib.request.urlopen(template).read().decode('UTF-8'), 'template')
    cmd.split_states('template')
    cmd.delete('template')
    cmd.split_chains(f'template_0001 and polymer.protein')
    cmd.delete(f'template_0001')    
    df = match_fasta_position(pdb2fasta(subject)[0], pdb2fasta(template))
    objs = list(cmd.get_object_list('all'))
    cmd.load(subject, 'subject')
    
    n=65  ##letter A
    for mol in objs:
        chain = cmd.get_chains(mol)[0]
        if chain in set(df['sequence'].tolist()):
            cmd.create(f'sub_{mol}', f"resi {int(df[df['sequence']==chain]['query_num'].min())}-{int(df[df['sequence']==chain]['query_num'].max())} and subject")
            cmd.alter(f'(sub_{mol})', f"chain='{chr(n)}'")
            cmd.fit(f'sub_{mol}', mol)
            n+=1
    cmd.save(f'{ac}.pdb', 'sub_*')

def build_homoligo_manual(subjects, template, ac='structure'):
    cmd.reinitialize()
    try:
        cmd.load(template, 'template')
    except:
        cmd.read_pdbstr(urllib.request.urlopen(template).read().decode('UTF-8'), 'template')
    cmd.split_states('template')
    cmd.delete('template')
    cmd.split_chains(f'template_0001 and polymer.protein')
    cmd.delete(f'template_0001')    
    objs = list(cmd.get_object_list('all'))
    
    n=65  ##letter A
    for chain, sub in subjects.items():
        df = match_fasta_position(pdb2fasta(sub)[0], [s for s in pdb2fasta('https://files.rcsb.org/view/7CQK.pdb') if s.id==chain])
        try:
            cmd.load(sub, 'subject')
        except:
            cmd.read_pdbstr(urllib.request.urlopen(sub).read().decode('UTF-8'), 'subject')
        mainChain = cmd.get_chains('subject')[0]
        cmd.create(f'sub_{chain}', f"resi {int(df[df['sequence']==chain]['query_num'].min())}-{int(df[df['sequence']==chain]['query_num'].max())} and subject and chain {mainChain}")
        cmd.alter(f'(sub_{chain})', f"chain='{chr(n)}'")
        cmd.fit(f'sub_{chain}', f'template_0001_{chain}')
        n+=1
        cmd.delete(f'subject')
    cmd.save(f'{ac}.pdb', 'sub_*')
    
def best_structure(ac):
    try:
        df = getfeatures(ac)
        leng = df['sequence.length'].iloc[0]*0.7
        pdb = pd.json_normalize([x  for x in df.dbReferences[0] if x['type']=='PDB'])
        pdb = pdb[pdb['properties.chains'].apply(lambda x: abs(eval(x.split('=')[-1]))) > leng]
        multi = pdb[pdb['properties.chains'].str.contains('/')]
        pdb = multi if not multi.empty else pdb
        pdb['properties.resolution'] = pdb['properties.resolution'].apply(lambda x: x.split(' ')[0]).astype(float)
        return pdb.loc[pdb['properties.resolution']==pdb['properties.resolution'].min()]['id'].iloc[0]
    except:
        return None
    
def getfeatures(uniprot_id):
    """return the Uniprot complete dataset for a given Uniprot ID"""
    try:
        r = requests.get("https://www.ebi.ac.uk/proteins/api/proteins?size=-1&accession=" + uniprot_id, 
                         headers={ "Accept" : "application/json"})
        data = pd.json_normalize(r.json())
        return data
    except:
        return str(uniprot_id) + ' not_found'


def planeEq(p1,p2,p3):
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    v1 = p3 - p1
    v2 = p2 - p1

    cp = np.cross(v2, v1)
    a, b, c = cp

    d = np.dot(cp, p3)
    return [a,b,c,d]

def lineEq(p1,p2):
    p1 = np.array(p1)
    p2 = np.array(p2)

    a, b, c = p2-p1

    return [a,b,c]

def angle_plane_line(a1,b1,c1,a2,b2,c2):
    return (a1*a2 + b1*b2 + c1*c2)/(sqrt(a1**2 + b1**2 + c1**2)*sqrt(a2**2 + b2**2 + c2**2))

def calc_run(run_file, dihedral):
    reactions = ['X1','X2', 'X3']
    CFCsBool = ['CFC_1_bool', 'CFC_2_bool', 'CFC_3_bool']
    pdbqt = os.path.basename(run_file[:-10])
    pdbqt = f"{os.path.dirname(os.path.dirname(run_file))}/{os.path.basename(run_file[:-10])}"

    run = os.path.basename(run_file)
    DIR = os.path.basename(os.path.dirname(run_file))
    ## PYMOL
    # try:
    cmd.reinitialize()
    cmd.load(run_file)
    # cmd.load(f"{pdbqt}.pdbqt")
    cmd.h_add(f'name {str(dihedral[2])} and not resn lys')
    n = cmd.count_atoms('not resn lys')

    angles = {}
    plane = cmd.get_coords(f"name {'+'.join(dihedral[1])} and not resn lys")
    a2,b2,c2,d2 = planeEq(*plane)
    
    # betalyase, decarboxylase
    for k,v in dihedral[0].items():
        line = cmd.get_coords(f'name {str(dihedral[2])}+{v} and not resn lys')
        a1,b1,c1 = lineEq(*line)
        angle = abs(angle_plane_line(a1,b1,c1,a2,b2,c2))
        angles.update({k:angle})
        
    if len(angles)==1:
        line = cmd.get_coords(f'name {str(dihedral[2])} or id {str(n-1)} and not resn lys')
        a1,b1,c1 = lineEq(*line)
        angle = abs(angle_plane_line(a1,b1,c1,a2,b2,c2))
        angles.update({'X1':angle})
        
    # aminotransferase
    line = cmd.get_coords(f'name {str(dihedral[2])} or id {str(n)} and not resn lys')
    a1,b1,c1 = lineEq(*line)
    angle = abs(angle_plane_line(a1,b1,c1,a2,b2,c2))
    angles.update({'X3':angle})
    
    d = cmd.get_distance(f'name {str(dihedral[3])} and not resn lys', 'name NZ and resn lys')
    results = pd.DataFrame([[DIR, os.path.basename(run_file), d, 
                            #  h_bond, 
                              *list(angles.values())]],
                            columns = ['DIR', 'run_file', 'd', 
                                      # 'h_bond',
                                       *list(angles.keys())]
                          )    
    results[CFCsBool] = results[['d', *reactions]].apply(lambda x: x[reactions]==max(x[reactions]), axis=1) 
    return results

def DLGdf(summary, runs_tab, resn):
    reactions = ['X1','X2', 'X3']
    CFCs = ['CFC_1', 'CFC_2', 'CFC_3']
    result = []
    pdbqt = summary.split('_PLP_')[1].split('_sum')[0]
    DIR = '_'.join(summary.split('_')[:-1])
    trg = list(enumerate(open(summary).readlines()))
    trg_file = [trg[n][1].split()[-1] for n, x in trg if 'Unpa' in x][0]

    conformations = []

    for idx, line in enumerate(open(summary).readlines()):
        if line.startswith('cluster'):
            clu = int(line.split(' ')[-1].strip())
        elif line.startswith('  runs'):
            runs = line.split()[1:]
        elif line.startswith('  energies'):
            energies = line.split()[1:]
            for en, run in zip(energies,runs):
                conformations.append([clu,int(run),float(en)])
    conformations = pd.DataFrame(conformations, columns=['rank', 'run', 'energy'])
    conformations['rank'] = conformations['rank'].astype(int)

    table_idx = [idx for idx, line in enumerate(open(summary).readlines()) 
                if line.endswith('+\n') 
                or line.startswith('Clustering information')
            ]

    df = pd.read_csv(summary, skiprows = table_idx[0]+1, nrows=table_idx[1]-table_idx[0]-3, 
                delim_whitespace=True, 
                names = ['rank', 'affinity', 'clust_rmsd','ref_rms', 'runs', 'rmsd_stdv','energy_stdv','best_run'],
                        # index_col = 'mode'
                )
    df = df.merge(conformations, on='rank')


    runs_tab['run'] = runs_tab['run_file'].apply(lambda x: int(x.replace('.pdbqt','')[-4:]))
    df = runs_tab.merge(df, on='run')
    clus = dict(df[['rank', 'run']].groupby('rank')['run'].apply(list))

    LC_runs = clus.get(df.sort_values(['runs','rank'], ascending=[False,True]).iloc[0]['rank'])
    BC_runs = clus.get(1)

    for k, r in clus.items():
        k_df = df[df['run'].isin(r)]
        result.append([pdbqt, 
                        len(k_df[k_df['CFC_1_bool']]), len(k_df[k_df['CFC_2_bool']]), len(k_df[k_df['CFC_3_bool']]), 
                        k_df['d'].astype(float).mean(), 
                        k_df['X1'].astype(float).mean(), k_df['X2'].astype(float).mean(), k_df['X3'].astype(float).mean(),
                        len(LC_runs), r,
                        ])

    return pd.DataFrame(result, 
                 columns = ['file', 
                            *CFCs,
                            'd mean',
                            f'X1_mean', f'X2_mean', f'X3_mean',
                            'LCC', 'run',
                           ]).explode('run').merge(df, on='run')
        

def plotCatalyitic(histo):
    CFCs = ['CFC_1', 'CFC_2', 'CFC_3']
    reactions = ['X1','X2', 'X3']
    
    plt.rcParams["hatch.linewidth"] = 8
    fig,ax = plt.subplots(figsize=(10,6))
    names = histo['file'].apply(lambda x: x.split('_')).tolist()[:2]
    titleName = '_'.join(sorted(set(names[0])&set(names[1]), key=names[0].index))

    plt.title(f"{titleName} catalytic plot", fontsize='x-large')

    histo = histo[[*CFCs, 'rank', 'runs', 'affinity']]
    histo['cat_runs'] = histo[CFCs].apply(list, axis=1)
    histo['reaction'] = [['decarboxylase', 'betalyase', 'aminotransferase']]*len(histo)
    histo = histo.explode(['reaction', 'cat_runs'])
    energies = histo.drop_duplicates('rank').sort_values('runs',ascending=True)

    mainRect = []
    colorsReaction = ['#467099ff', '#4e9976ff', '#99b802ff']
    off = 0.15

    for n, (i, en) in enumerate(energies.iterrows()):
        if all([en['affinity'] == energies['affinity'].min(), # Best cluster and largest cluster
                en['runs'] == energies['runs'].max()]):
            edgecolor = (.9, .8, .1)
            facecolor = (.8, .1, .1)
            alpha=.8
            # facecolor = pd.DataFrame([[.55, .8, .9], [1, .4, .25]]).mean().tolist()
            hatch = r"/"
            # hatch = r""
        elif en['affinity'] == energies['affinity'].min(): #  Best cluster
            edgecolor = 'black'
            facecolor = (.8, .1, .1)
            alpha=.8
            hatch = ''
        elif en['runs'] == energies['runs'].max(): #  Largest cluster
            edgecolor = 'black'
            facecolor = (.9, .8, .1)
            alpha=.8
            hatch = ''
        else:
            edgecolor = 'black'
            facecolor = (.9, .9, .9)
            alpha=1
            hatch = ''

        rect = patches.Rectangle((en['affinity']-0.16, 0), 0.32, en['runs'], 
                                linewidth=1,
                                edgecolor=edgecolor, 
                                facecolor=facecolor,
                                 hatch=hatch
                                )
        ax.add_patch(rect)
        off = 0.15
        for r, color in zip(CFCs, colorsReaction):
            rect = patches.Rectangle((en['affinity']-off, 0), 0.1, en[r], 
                                    linewidth=0.5,
                                    edgecolor='gray', 
                                    facecolor=color
                                    )
            ax.add_patch(rect)
            off-=0.1
        
    ax.set_xlim(energies['affinity'].min()-0.5,energies['affinity'].max()+0.5)
    ax.set_ylim(0,energies['runs'].max()+10)
    ax.set_xlabel('Binding Energy (kcal/mol)', fontsize= 'large')
    ax.set_ylabel('Runs', fontsize= 'large')
    ax.legend(['all', 'Cα-COOH', 'Cα-Cβ', 'Cα-Hα'],
              fancybox=True, 
              framealpha=1, 
              shadow=True, 
              borderpad=1
             )
    plt.close()
    return fig, ax

def calc_ade(run, trg):
    cwd = os.getcwd()
    lig = pd.DataFrame(map(lambda x: x.split(),re.findall('ATOM.*', open(run).read()))) # read only HETATM lines
    atoms = lig[~lig[3].isin(['PLP', 'LYS'])][1].astype(int).tolist() # atom without PLP and lys
    os.chdir(os.path.dirname(run))
    ade = os.popen(f'~/G_BIOSCIENZE/ADFRsuite_x86_64Linux_1.0/bin/pythonsh ~/G_BIOSCIENZE/ADFRsuite_x86_64Linux_1.0/CCSBpckgs/ADFR/bin/ade.py -t ../{trg} -l {os.path.basename(run)}').read()
    os.chdir(cwd)

    idx = [StringIO(ade[x.span()[0]:x.span()[1]]) for x in list(re.finditer('_+.*\n(\s+\d+.*\n)+', ade))]
    lrr, ll = [pd.read_csv(x, skiprows = 1, header=None,keep_default_na=False,
                    delim_whitespace=True) for x in [idx[0], idx[2]]]
    lrr[0] = lrr[0]+1  ### numero atomo sbagliato in LRR, la numerazione parte da 0 anziché da 1
    lrr_energy = lrr[lrr[0].isin(atoms)][2].sum()
    ll[[1,10]] = ll[1].str.split('-',expand=True).astype(int)
    ll_energy = ll[(ll[1].isin(atoms))|(ll[10].isin(atoms))][5].sum()
    energy = ll_energy + lrr_energy
    return [os.path.basename(run), lrr_energy]