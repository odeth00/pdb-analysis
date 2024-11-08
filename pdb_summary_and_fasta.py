#!/home/odeth/miniforge/bin/python
#!/home/alejandro/bin/yopython
#!/home/alejandro/miniconda3/envs/yo/bin/python
# Your shebang at the beginning

# ====================== imports ====================

import sys
import pandas as pd
import warnings
import argparse
from io import StringIO
import textwrap

warnings.simplefilter(action='ignore', category=FutureWarning)

# =================== functions ====================

amino_acids = [
    "ALA", # Alanine
    "ARG", # Arginine
    "ASN", # Asparagine
    "ASP", # Aspartic acid
    "CYS", # Cysteine
    "GLN", # Glutamine
    "GLU", # Glutamic acid
    "GLY", # Glycine
    "HIS", # Histidines
    "ILE", # Isoleucine
    "LEU", # Leucine
    "LYS", # Lysine
    "MET", # Methionine
    "PHE", # Phenylalanine
    "PRO", # Proline
    "SER", # Serine
    "THR", # Threonine
    "TRP", # Tryptophan
    "TYR", # Tyrosine
    "VAL"  # Valine
]

one_to_three = {
    'A': 'ALA',
    'R': 'ARG',
    'N': 'ASN',
    'D': 'ASP',
    'C': 'CYS',
    'E': 'GLU',
    'Q': 'GLN',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'L': 'LEU',
    'K': 'LYS',
    'M': 'MET',
    'F': 'PHE',
    'P': 'PRO',
    'S': 'SER',
    'T': 'THR',
    'W': 'TRP',
    'Y': 'TYR',
    'V': 'VAL'
}

three_to_one = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',
    '---':'-'
}

def pdb_atoms2df(pdb_path):
    """
    Reads a pdb file from pdb_path, and returns a dataframe with all the ATOM, HELIX, and SHEET lines.
    """
    atom_lines = []
    helix_lines = []
    sheet_lines = []

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atom_lines.append(line)
            elif line.startswith('HELIX'):
                helix_lines.append(line)
            elif line.startswith('SHEET'):
                sheet_lines.append(line)
    
    atom_df = parse_atoms(atom_lines)
    
    helix_df = parse_helix(helix_lines)
    
    sheet_df = parse_sheet(sheet_lines)

    return atom_df, helix_df, sheet_df

def parse_atoms(atom_lines):
    """
    Parse atom lines in a pdb file
    """
    atom_lines_text = ''.join(atom_lines)
    atom_lines_io = StringIO(atom_lines_text)

    col_specs = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), 
                (22, 26), (26, 27), (30, 38), (38, 46), (46, 54), (54, 60), 
                (60, 66), (76, 78), (78, 80)]
    col_names = ['record_name', 'serial', 'name', 'altLoc', 'resName', 'chainID', 
                'resSeq', 'iCode', 'x', 'y', 'z', 'occupancy', 'tempFactor', 
                'element', 'charge']

    df = pd.read_fwf(atom_lines_io, colspecs=col_specs, names=col_names, header=None)

    numeric_columns = ['serial', 'resSeq', 'x', 'y', 'z', 'occupancy', 'tempFactor']
    df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors='coerce')

    return df

def parse_helix(helix_lines):
    """
    Parse helix lines in a pdb file
    """
    helix_lines_text = ''.join(helix_lines)
    helix_lines_io = StringIO(helix_lines_text)

    col_specs = [(7, 10), (11, 14), (15, 18), (19, 20), 
                 (21, 25), (25, 26), (27, 30), (31, 32), 
                 (33, 37), (37, 38), (38, 40), (40, 70), 
                 (71, 76)]
    col_names = ['serNum', 'helixID', 'initResName', 'initChainID', 
                 'initSeqNum', 'initICode', 'endResName', 'endChainID', 
                 'endSeqNum', 'endICode', 'helixClass', 'comment', 
                 'length']

    helix_df = pd.read_fwf(helix_lines_io, colspecs=col_specs, names=col_names, header=None)

    numeric_columns = ['serNum', 'initSeqNum', 'endSeqNum', 'helixClass', 'length']
    helix_df[numeric_columns] = helix_df[numeric_columns].apply(pd.to_numeric, errors='coerce')
    
    return helix_df

def parse_sheet(sheet_lines):
    """
    Parse sheet lines in a pdb file
    """
    sheet_lines_text = ''.join(sheet_lines)
    sheet_lines_io = StringIO(sheet_lines_text)

    col_specs = [(7, 10), (11, 14), (14, 16), (17, 20), (21, 22), 
                 (22, 26), (26, 27), (28, 31), (32, 33), (33, 37), 
                 (37, 38), (38, 40), (41, 45), (45, 48), (49, 50), 
                 (50, 54), (54, 55), (56, 60), (60, 63), (64, 65), 
                 (65, 69), (69, 70)]
    
    col_names = ['strand', 'sheetID', 'numStrands', 'initResName', 'initChainID', 
                 'initSeqNum', 'initICode', 'endResName', 'endChainID', 
                 'endSeqNum', 'endICode', 'sense', 'curAtom', 'curResName', 
                 'curChainId', 'curResSeq', 'curICode', 'prevAtom', 'prevResName', 
                 'prevChainId', 'prevResSeq', 'prevICode']

    sheet_df = pd.read_fwf(sheet_lines_io, colspecs=col_specs, names=col_names, header=None)

    numeric_columns = ['strand', 'numStrands', 'initSeqNum', 'endSeqNum', 'sense', 'curResSeq', 'prevResSeq']
    sheet_df[numeric_columns] = sheet_df[numeric_columns].apply(pd.to_numeric, errors='coerce')
    
    return sheet_df 

def get_final_composition_4(df):
    df4 = df.loc[(df['name'] == 'CA'),['chainID','resName'] ]
    df5 = pd.crosstab(index=df4['resName'], columns=df4['chainID'])
    return df5

    
def mean_tempf(df):
    """
    Obtain the Bfactor mean group by chaiID
    """
    df1 = df.groupby('chainID')['tempFactor'].mean()
    return df1    


def missing_resSeq(df):
    """
    Finds missing resSeq numbers in each chain.
    returns a df.Series where chains are index and values are lists
    """
    chains = df['chainID'].unique()
    missing = {}
    for c in chains:
        missing[c] = []
        nums = df.loc[(df.chainID==c) & (df.name == 'CA'), 'resSeq'].tolist()
        for num in range(nums[0], nums[-1]+1):
            if not num in nums:
                missing[c].append(str(num))
    return pd.Series(missing)

def start_and_end(df):
    """
    Obtain start and end of resSeq group by chainID
    """
    chain_start = df.groupby('chainID')['resSeq'].agg(['min'])
    chain_end = df.groupby('chainID')['resSeq'].agg(['max'])
    ## yo
    print(type(chain_start), chain_start)
    return chain_start, chain_end
    
def chain_lens(df): 
    """
    Obtain length of every chain and add columns with the temp factor and start and end resSeq also missing resSeq
    """
    ca = df[df['name'] == 'CA']
    ca_chains = ca.groupby('chainID').size() 
    start_end_chain = start_and_end(df)
    ca_chains_df = ca_chains.to_frame(name='Residue_Count')
    chain_start, chain_end = start_end_chain 
    ca_chains_df['Start'] = chain_start
    ca_chains_df['End'] = chain_end
    temp_mean = mean_tempf(df) 
    ca_chains_df['TempF_Mean'] = temp_mean
    ca_chains_df['missing_resSeq'] = missing_resSeq(df)
    return ca_chains_df

def fasta_format(df,pdb):
    """
    Obtain fasta format of every chain
    """
    pdb_name = pdb.split('/')[-1][0:-4]
    secuencias_aa = {}
    for chain in df['chainID'].unique():
        df_cadena = df.loc[df['chainID'] == chain]
        resNames = df_cadena.loc[df_cadena['name'] == 'CA',['resSeq','resName']]
        
        complete_resSeq = pd.Series(range(resNames['resSeq'].min(), resNames['resSeq'].max() + 1), name='resSeq')
        df_complete = pd.merge(complete_resSeq, resNames, on='resSeq', how='left')
        df_complete['resName'].fillna('---', inplace=True)
        
        secuencia_aa = ''.join(three_to_one[res] for res in df_complete['resName'])
        secuencias_aa[chain] = secuencia_aa
        
    for chain, sequence in secuencias_aa.items():
        sequence = textwrap.fill(sequence, width=80, break_on_hyphens=False)
        print(f">{pdb_name}{chain}\n{sequence}")
    
    return secuencias_aa

def do_one(atom_df, helix_df, sheet_df, pdb):
    chlen = chain_lens(atom_df) 
    print("===== Chain summary ====")
    print(chlen)

    comp = get_final_composition_4(atom_df) 
    print("===== AA composition ====")
    print(comp)

    print("===== Fasta sequences ====")
    fasta_format(atom_df, pdb) 

# ======================== main =========================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Makes a summary of a PDB file '.pdb'")
    parser.add_argument('-f', '--fasta_only', action='store_true', help='Enables fasta only mode.')
    parser.add_argument('pdb', help='File path')
    args = parser.parse_args()
    
    pdb = args.pdb
    atom_df, helix_df, sheet_df = pdb_atoms2df(pdb)

    if args.fasta_only:
        fasta_format(atom_df, pdb)
        sys.exit()

    do_one(atom_df, helix_df, sheet_df, pdb)

# ============ end =============
