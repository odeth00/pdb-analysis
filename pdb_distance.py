#!/home/odeth/miniforge/bin/python
#!/home/alejandro/bin/yopython
#!/home/odeth/miniforge/envs/jupyter/bin/python
# Your shebang at the beginning

# %% ======================== imports ===========================
import sys
from glob import glob
import pandas as pd
import argparse
from io import StringIO
from scipy.spatial import distance
from scipy.spatial.distance import cdist
from collections import defaultdict

# %% ======================== data ===============================
# Define the mapping from three-letter to one-letter amino acid codes
amino_acid_mapping = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# Create a defaultdict with default value 'X' for missing keys
one_let = defaultdict(lambda: 'X', amino_acid_mapping)


# %% ========== functions ============

def pdb_atoms2df(pdb_path):
    """
    Reads a pdb file from pdb_path, and returns a df with all the ATOM lines.
    """
    atom_lines = []
    with open(pdb_path,'r') as f:
        for line in f:
            if line.startswith('ATOM'):
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
    df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors='coerce')

    return df

def get_user_chain_off(df):
    """
    Filter by type chain
    """
    chains = df['chainID'].unique()
    if len(chains) > 1:
        wanted_chains = input(f"Enter the chain IDs you're interested in (Available chains: {chains}): ").upper().strip()
        df = df[df['chainID'].isin(list(wanted_chains))]
    else:
        print(f"Only one chain available: {chains[0]}")
    return df


def reference_atom_off(df):
    ok_chains = df['chainID'].unique()
    reference_atoms = {}
    for chain in ok_chains:
        min_serial = df.loc[df.chainID == chain,'serial'].min()
        max_serial = df.loc[df.chainID == chain,'serial'].max()
        wanted_atoms = input(f"Enter the serial numbers of the reference atoms for chain {chain} ({min_serial}-{max_serial}): ").strip()
        reference_atoms[chain] = [int(atom) for atom in wanted_atoms.split() if atom.isdigit()]

    return reference_atoms

def get_reference(df):
    """
    Para cada cadena, deja al usuario escoger ciertos residuos
    Convierte los residuos a atomos.
    Regresa
    refrenece_atoms : los atomos de esos residuos
    reference_residues : un dict cadena: str("los residuos elegidos")
    """
    ok_chains = df['chainID'].unique()
    reference_atoms = set()
    reference_residues = {}
    for chain in ok_chains:
        min_serial = df.loc[df.chainID == chain,'resSeq'].min()
        max_serial = df.loc[df.chainID == chain,'resSeq'].max()
        answer = input(f"Enter the numbers of the reference residues for chain {chain} ({min_serial}-{max_serial}) or 'enter' for None: ").strip()
        reference_residues[chain]=answer
        if answer == "":
            continue
        wanted = set()
        parts = answer.split()
        for part in parts:
            if '-' in part:
                start,end  = [int(x) for x in part.split('-')]
                wanted.update(list(range(start,end+1)))
            else:
                wanted.add(int(part))
        for resSeq in wanted:
            atoms = df.loc[(df.chainID == chain) & (df.resSeq == int(resSeq)), 'serial'].tolist()
            reference_atoms.update(atoms)

    return reference_atoms, reference_residues

def calculate_distances_off(df, selected_atoms, max_distance):
    atoms_df = df[df['serial'].isin(selected_atoms)]
    for i, atom1 in atoms_df.iterrows():
        for j, atom2 in df.iterrows():
            if i < j:
                dist = distance.euclidean(atom1[['x', 'y', 'z']], atom2[['x', 'y', 'z']])
                if dist <= max_distance:
                    print(f"Distance between atom {atom1['serial']} and {atom2['serial']}: {dist:.2f} Å")

def get_close_atoms(df, wanted_atoms, max_dist):
    """
    Para wanted_atoms calcula sus distancias a todos los demas atomos
    Y luego escoge los que tienen distancia menor a max_dist
    Regresa:
    pairs (list of tuples): cada tupla tiene (referece_atom, a_near_atom)
    """
    df1 = df.loc[df.serial.isin(wanted_atoms)] # rebanada de df cuyos rows matchean con wanted_atoms
    df2 = df  # Por claridad
    # Extract coordinates
    coords1 = df1[['x', 'y', 'z']].values
    coords2 = df2[['x', 'y', 'z']].values

    # Usamos cdist en vez de distance, pues cdist hace un set vs otro.
    distances = cdist(coords1, coords2, metric='euclidean')

    # Convert to DataFrame
    distance_df = pd.DataFrame(dis