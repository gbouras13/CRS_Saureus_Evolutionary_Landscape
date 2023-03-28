"""
The snakefile that runs the pipeline.
Manual launch example:


"""

import os
import pandas as pd

### DEFAULT CONFIG FILE

BigJobMem = 4000
SmallJobMem = 100
BigJobCpu = 2

Plasmid_csv = 'ghais_plasmids.csv'
BAKTA_DB = '/Volumes/VERBATIM_HD_1/bakta_db/db'

### DIRECTORIES

include: "rules/directories.smk"

# get if needed
INPUT = config['Input']
OUTPUT = config['Output']

## parse file from csv of 
PLASMIDS = pd.read_csv(Plasmid_csv, header=None).loc[:, 0].tolist()
print(PLASMIDS)

# Import rules and functions
include: "rules/targets_mash.smk"
include: "rules/mash.smk"
include: "rules/bakta.smk"
include: "rules/panaroo.smk"
include: "rules/jaccard.smk"
include: "rules/abricate.smk"

rule all:
    input:
        TargetFilesMash