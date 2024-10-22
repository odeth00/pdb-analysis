# pdb-analysis
This suite contains the following programs to analyze pdb files.
Send comments and bugs to Odeth Valencia.
Some scripts allow the user to interactively define what is wanted.

## pdb_to_tabular.py (non-interactive)
USAGE: pdb_to_tabular.py [-h] [-s] pdb outfile
This script converts PDB files to TSV or XLSX format.
By default, it only reports the ATOM section. With option '-s' it annotates HELIX or SHEET in an additional column.

## pdb_summary_and_fasta.py (non-interactive)
UASGE: pdb_summary_and_fasta.py [-h] [-f] pdb
NOTE: Prints to STDOUT
This script summarizes PDB files. 
For each chain, it reports chain_name, number_of_residues, first and last residue, mean bfactor, 
    list of missing residues, amino acid composition and fasta sequence
With -f, it only reports the fasta sequences.

## pdb_filter.py (interactive)
USAGE: pdb_filter.py [-h] [-c] [-a] [-r] [-n] pdb outfile.pdb
NOTE1: All options can be combined. Example: -can (filter by chain, atom and residue number)
NOTE2: If none of -[carn] is given, runs without interaction and extracts all ATOM lines.
This script provides a command-line interface for filtering the ATOM section PDB files based 
on criteria such as amino acid chain, residue type, atom type and residue number.
It offers interactive filtering options, and saves the filtered results back into PDB format.

## pdb_distance.py (interactive)
USAGE: pdb_distance.py [-h] pdb_path outfile
This script analyzes distances between atoms in a PDB file. 
It reads the PDB file, prompts the user to select reference residues, and specify a maximum distance. 
It then calculates distances between corresponding atoms and any other atom in the pdb.
Finally, it prints atom pairs with distances below the specified threshold.

## pdb_filter_trajectory.py (interactive)
USAGE: pdb_filter_trajectory.py [-h] pdb outdir
This script is designed to extract individual PDB files from a multi-PDB trajectory file.
It parses the trajectory file, identifies each model, and prompts the user to specify which models they want to extract.
The user can input specific model numbers, ranges of models, or define a step size for jumping through models.
The script then saves the selected models as separate PDB files in the specified output directory.
