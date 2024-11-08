#!/home/odeth/miniforge/bin/python
#!/home/alejandro/bin/yopython
#!/home/alejandro/mambaforge/envs/yo/bin/python
# Your shebang at the beginning

# %% ================= imports =================

import sys
import os
import re
from glob import glob
import argparse
"""
# Para acordarnos
here = os.getcwd()
sys.path.insert(0, here)
# import resumen_pdb as rpdb
from resumen_pdb import *
"""

# ==================== functions ===============

def parse_wanted_files(wanted, max):
    wanted = wanted.strip().split()
    numbers = []
    for w in wanted:
        try:
            w = int(w)
            numbers.append(w)
            continue
        except:
            pass
        if match:=re.match('(\d+)-(\d+)',w):
            start = int(match.group(1))
            end = int(match.group(2))
            for i in range(start,end+1):
                numbers.append(i)
        elif  match:=re.match('j(\d+)',w):
            step= int(match.group(1))
            for i in range(0,max,step):
                numbers.append(i)
        else:
            print(f"Parsing {w} not implemented")
    return numbers

# ========== main =========
    
if __name__ == '__main__':

    epi="""
    Separate filenums, ranges, and jumps with spaces:
    ranges are 'start-end' (example: 10-18);
    jumps are 'j[jump_size]' (example: j5)
    """
