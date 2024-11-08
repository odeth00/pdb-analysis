#!/home/odeth/miniforge/bin/python
#!/home/alejandro/bin/yopython
#!/home/alejandro/miniconda3/envs/yo/bin/python
# Your shebang at the beginning


# ====================== imports ====================

import pandas as pd
import numpy as np
import re
import os
import sys
import time
from collections import defaultdict
from io import StringIO
import argparse

# ==================== functions ==================

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
    
def get_atoms(pdb):
    """
    Search the ATOM lines of a PDB, extract the fields, and add them as new rows to a pd.DataFrame. 
    Return the DataFrame
    """
    start = time.time()
    df = pd.DataFrame()
    with open(pdb,'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                campos = parse_pdb_atom_line(line)
                
                new_row = pd.DataFrame([campos])
    
                df = pd.concat([df, new_row], ignore_index=True)

    end = time.time()
    print(f"get_atoms took {end-start} segundos")
    
    return df


def parse_pdb_atom_line(line):
    """ 
    Extract the fields from an 'ATOM' line of the PDB.
    Return a dictionary field:value.    
    """
    parsed_data = {
        'record_name': line[0:6].strip(),
        'serial': int(line[6:11].strip()),
        'name': line[12:16].strip(),
        'altLoc': line[16].strip(),
        'resName': line[17:20].strip(),
        'chainID': line[21].strip(),
        'resSeq': int(line[22:26].strip()),
        'iCode': line[26].strip(),
        'x': float(line[30:38].strip()),
        'y': float(line[38:46].strip()),
        'z': float(line[46:54].strip()),
        'occupancy': float(line[54:60].strip()),
        'tempFactor': float(line[60:66].strip()),
        'element': line[76:78].strip(),
        'charge': line[78:80].strip()
    }
    return parsed_data

# ======================= main =======================

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description= "Convert a pdb to tabular format'.pdb'" )
    parser.add_argument('pdb', help='File path')
    parser.add_argument('outfile', help="Output file ('.tsv' or '.xlsx)")
    parser.add_argument('-s', action='store_true', help= 'Add column with secondary structure')

    args = parser.parse_args()

    pdb = args.pdb
    
    atom_df, helix_df, sheet_df = pdb_atoms2df(pdb)
    if args.s:
        atom_df['SStructure'] = ''
        for _, row in helix_df.iterrows():
            mask = (
                (atom_df['chainID'] == row['initChainID']) &
                (atom_df['resSeq'] >= row['initSeqNum']) &
                (atom_df['resSeq'] <= row['endSeqNum'])
            )
            atom_df.loc[mask, 'SStructure'] = 'H'
        
        for _, row in sheet_df.iterrows():
            mask = (
                (atom_df['chainID'] == row['initChainID']) &
                (atom_df['resSeq'] >= row['initSeqNum']) &
                (atom_df['resSeq'] <= row['endSeqNum'])
            )
            atom_df.loc[mask, 'SStructure'] = 'B'

    df = atom_df

 if args.outfile.endswith('xlsx'):
        with pd.ExcelWriter(f'{args.outfile}') as writer: 
            for chain_id in df['chainID'].unique():
                chain_df = df[df['chainID'] == chain_id]
                chain_df.to_excel(writer, sheet_name=str(chain_id), index=False)
        print(f"Saved ATOM data {args.outfile}")

    elif args.outfile.endswith('tsv'):
        atom_df.to_csv(args.outfile, sep='\t', index=False)
        print(f"Saved ATOM data as {args.outfile}")

    else:
        print(f"I dont know the format for {args.outfile}")

    sys.exit()
