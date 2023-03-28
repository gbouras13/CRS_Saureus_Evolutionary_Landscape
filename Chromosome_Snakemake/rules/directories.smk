"""
Database and output locations 
Ensures consistent variable names and file locations for the pipeline.
"""

import os


### OUTPUT DIRECTORY
if config['Output'] is None:
    OUTPUT = 'Assemblies_Output'
else:
    OUTPUT = config['Output']

# REFERENCE_DIR
REFERENCE = "Reference"
# pharokka db
PHAROKKA_DB="/hpcfs/users/a1667917/pharokka_db"



### OUTPUT DIRs
FLAGS = os.path.join(OUTPUT, 'FLAGS')




# ghais dirs
SNIPPY_PAIR = os.path.join(OUTPUT, 'SNIPPY_PAIR')
NUCDIFF = os.path.join(OUTPUT, 'NUCDIFF')
MLST = os.path.join(OUTPUT, 'MLST')
PANAROO = os.path.join(OUTPUT, 'PANAROO')
IQTREE = os.path.join(OUTPUT, 'IQTREE')
PHISPY = os.path.join(OUTPUT, 'PHISPY')
ISESCAN = os.path.join(OUTPUT, 'ISESCAN')
ABRICATE = os.path.join(OUTPUT, 'ABRICATE')
SNIFFLES = os.path.join(OUTPUT, 'SNIFFLES')

# needs to be created prior
if not os.path.exists(SNIPPY_PAIR):
  os.makedirs(SNIPPY_PAIR)

if not os.path.exists(NUCDIFF):
  os.makedirs(NUCDIFF)








