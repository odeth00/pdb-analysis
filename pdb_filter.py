#!/home/odeth/miniforge/bin/python
#!/home/alejandro/bin/yopython
#!/home/alejandro/Dropbox/work/parse_pdb/pdb_env/bin/python
#!/home/alejandro/miniconda3/envs/yo/bin/python
# Your shebang at the beginning

import sys
import os
import re
import pandas as pd
import argparse
import subprocess
from io import StringIO

uno = {
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
    'VAL': 'V'
}

tres = {
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

# ===== functions ========

def rebuild_pdb(filtered, original):
    kept_atoms = set()
    rebuilt = filtered.replace('.pdb','_rebuilt.pdb')
    with open(filtered, 'r') as fh:
        for line in fh:
            line = line.replace(' ','')
            kept_atoms.add(line)
    with open(original,'r') as ori, open(rebuilt,'w') as out:
        for line in ori:
            if not line.startswith('ATOM'):
                out.write(line)
            elif line.replace(' ','') in kept_atoms:
                out.write(line)
    print(f"Wrote {rebuilt}", file=sys.stderr)

def pdb_atoms2df(pdb_path):
    """ 
    Reads a pdb file from pdb_path, and returns a df with all the ATOM lines.
    """
    # Read atom lines
    atom_lines = []
    with open(pdb_path,'r') as f:
        for line in f:
            if line. startswith('ATOM'):
                atom_lines.append(line)
    
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
    df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors='coerce').fillna(0)

    return df

def write_pdb_from_df(df, filename):
    """
    Create a PDB from a DataFrame
    """
    # print(df)
    df = df.fillna('')  # Replace NaN with empty strings
    # print(df.dtypes) # To debug wrong types
    with open(filename, 'w') as f:
        for index, row in df.iterrows():
            f.write(
                    f"{row['record_name']:<6}"
                    f"{row['serial']:5}  "
                    f"{row['name']:3}"
                    f"{row['altLoc']:1}"
                    f"{row['resName']:3} "
                    f"{row['chainID']:1}"
                    f"{row['resSeq']:4}"
                    f"{row['iCode']:1}   "
                    f"{row['x']:8.3f}"
                    f"{row['y']:8.3f}"
                    f"{row['z']:8.3f}"
                    f"{row['occupancy']:6.2f}"
                    f"{row['tempFactor']:6.2f}          "
                    f"{row['element']:>2}"
                    f"{row['charge']:2}\n"
                    )


def read_atoms_pdb(pdb):
    
    if pdb.endswith('.tsv'):
        file_p = pdb
    else:
        file_p = pdb + '.tsv'
    print("Trying to open file:", file_p)
    if not os.path.exists(file_p):
        subprocess.run(['python', 'pdb_to_csv.py', pdb])
    df = pd.read_csv(file_p, sep='\t')
    df = df.drop(df.columns[0], axis=1)
    
    return df

options  = {'chainID' : ['A','B','C','D','E'], 'resName' : 'ACDEFGHIKLMNPQRSTVWY', 'name' : ['CA','CB', 'C', 'O','N'], 'resSeq': {'A': '4-207','B':'7-98', 'C':'12-300','D':'32-167', 'E':'8-200'}}

def get_user_chain(df):
    """
    Filter by type chain
    """
    valid = df['chainID'].unique()
    wanted = input((f"Choose 'chainID's, separated by spaces. 'Enter' to choose all. ({valid}) ")).upper().strip()
    
    if wanted == '*' or wanted == '':
        wanted = ''.join(valid)
        return df,wanted.split()

    wanted = [x for x in list(wanted) if x in valid]
    # filtro = f'df.loc[df.chainID.isin({wanted})]'
    # print(filtro)
    # df = eval(filtro)
    df = df.loc[df.chainID.isin(wanted)]
    return df, wanted

def get_user_resName(df):
    """
    Filter by specific residues or aminoacids 
    """
    valid =  df['resName'].unique()
    valid_uno = ''.join( sorted([uno[x] for x in valid]) )
    wanted = input((f"Choose 'resName's. 'Enter' to choose all. ({valid_uno}) ")).upper().strip()

    if wanted == '*' or wanted == '':
        wanted = valid
        return df
    wanted_list = [tres[x] for x in wanted if x in tres]
    # filtro = f'df.loc[df.resName.isin({wanted_list})]'
    # print(filtro)
    # df = eval(filtro)
    df = df.loc[df.resName.isin(wanted_list)]
    return df

def get_user_name(df):
    """
    Filter by kind of atoms
    """
    valid = sorted(df['name'].unique())
    wanted = input((f"Choose atoms types, separated by spaces. 'Enter' to choose all types: ({valid}) ")).upper().strip()
    if wanted == '*' or wanted == '':
        wanted = valid
    else:
        wanted = wanted.split()
    # filtro = f'df.loc[df.name.isin({wanted})]'
    # print(filtro)
    # df = eval(filtro)
    df = df.loc[df.name.isin(wanted)]
    return df

def get_user_resSeq(df):
    """
    Filter by specific residue sequence
    """
    ok_chains = df['chainID'].unique()
    
    pattern = re.compile(r'^\d+(-\d+)?$')
    dfs = []
    for chain in ok_chains:
        min_res = df.loc[df.chainID == chain,'resSeq'].min()
        max_res = df.loc[df.chainID == chain,'resSeq'].max()
        answer = input(f"ResNums for chain {chain}. Separate with space or comma. Use '-' for ranges. 'Enter' to choose all. (Valid are: {min_res}-{max_res}) ").strip().lower()
        if answer == '' or answer == '*':
            # filtro = f"df.loc[df.chainID == '{chain}']"
            # tmp = eval(filtro)
            tmp = df.loc[df.chainID==chain]
            dfs.append(tmp)
            
        else:
            wanted_list = re.split('[\s,]+',answer)
            wrong = [x for x in wanted_list if not pattern.match(x)]
            if wrong:
                sys.exit(f"Invalid residue indicators : {wrong}")
            parsed_wanted_list = set()
            
            for val in wanted_list:
                if '-' in val:
                    start, end = [int(x.strip()) for x in val.split('-')]  # Strip and convert each part to int
                    rango = list(range(start, end + 1))
                    parsed_wanted_list.update(rango)  # Extend the new list with the range
                else:
                    parsed_wanted_list.add(int(val))  # Convert to int and append to the new list
            
            parsed_wanted_list = sorted(parsed_wanted_list)
            
            # filtro = f"df.loc[(df.chainID == '{chain}') & (df.resSeq.isin({parsed_wanted_list}))]"
            # tmp = eval(filtro)
            tmp = df.loc[(df.chainID == chain) & (df.resSeq.isin(parsed_wanted_list))]
            dfs.append(tmp)
    df = pd.concat(dfs)
    return df

# ========== main ============

if __name__ == '__main__':

    # Initialize the parser
    parser = argparse.ArgumentParser(description='Filters a pdb by several criteria')

    # Define arguments
    parser.add_argument('-c', '--chain', help='Filter by chain', action='store_const', const='chainID')
    parser.add_argument('-a', '--atom', help='Filter by type of atom (CA, N, O, etc.)', action='store_const', const='name')
    parser.add_argument('-r', '--residue', help='Filter by type of residue (ALA, ARG, etc)', action='store_const', const='resName')
    parser.add_argument('-n', '--number', help='Filter by residue number', action='store_const', const='resSeq')
    parser.add_argument('pdb', help='PDB file')
    parser.add_argument('outfile', help='Outfile.pdb')

    # Parse arguments
    args = parser.parse_args()

    # Convert args namespace to a dictionary
    args_dict = vars(args)

    # Fields wanted
    fields = [x for x in ['chain','atom','residue','number'] if args_dict[x]]



