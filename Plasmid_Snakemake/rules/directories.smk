"""
Database and output locations for Hecatomb
Ensures consistent variable names and file locations for the pipeline.
"""


### OUTPUT DIRECTORY
if config['Output'] is None:
    OUTPUT = 'Plasmid_Analysis_Output'
else:
    OUTPUT = config['Output']


### OUTPUT DIRs
LOGS = os.path.join(OUTPUT, 'LOGS')
RESULTS = os.path.join(OUTPUT, 'RESULTS')
PLASDB = os.path.join(RESULTS, 'PLASDB')
TMP = os.path.join(OUTPUT, 'TMP')
PLASMIDS_EXTRACTED = os.path.join(OUTPUT, 'PLASMIDS_EXTRACTED')
PANAROO = os.path.join(OUTPUT, 'PANAROO')
#PROKKA = os.path.join(TMP, 'PROKKA')
BAKTA = os.path.join(TMP, 'BAKTA')
MASH = os.path.join(OUTPUT, 'MASH')
MOBTYPER = os.path.join(OUTPUT, 'MOBTYPER')
PLASMID_GFFS = os.path.join(OUTPUT, 'PLASMID_GFFS')
SUMMARIES = os.path.join(OUTPUT, 'SUMMARIES')
ABRICATE =  os.path.join(OUTPUT, 'ABRICATE')


# for major conserved plasmid 

MAFFT = os.path.join(OUTPUT, 'MAFFT')
MAFFT_RENAME = os.path.join(MAFFT, 'MAFFT_RENAME')


#needs to be created before plasdb is run
if not os.path.exists(PLASDB):
    os.makedirs(PLASDB)






