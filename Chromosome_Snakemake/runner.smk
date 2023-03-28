"""
The snakefile that runs the pipeline.
Manual launch example:
"""

import os

### DEFAULT CONFIG FILE

configfile: os.path.join(  'config', 'config.yaml')


### DIRECTORIES

include: "rules/directories.smk"

# get if needed
CSV = config['csv']
OUTPUT = config['Output']
BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]
SmallJobMem = config["SmallJobMem"]


# Parse the samples and read files
include: "rules/samples.smk"
dictReads = parseGess(CSV)
SAMPLES = list(dictReads.keys())

# define functions 

def get_T0_r1(wildcards):
    return dictReads[wildcards.sample]["T1_R1"]
def get_T0_r2(wildcards):
    return dictReads[wildcards.sample]["T1_R2"]
def get_T0_gbk(wildcards):
    return dictReads[wildcards.sample]["T0_gbk"]
def get_T0_fasta(wildcards):
    return dictReads[wildcards.sample]["T0_fasta"]
def get_T1_fasta(wildcards):
    return dictReads[wildcards.sample]["T1_fasta"]
def get_T1_gbk(wildcards):
    return dictReads[wildcards.sample]["T1_gbk"]


def get_all_T0_gffs(wildcards):
    return [dictReads[s]["T0_gff"] for s in SAMPLES]
def get_all_T1_gffs(wildcards):
    return [dictReads[s]["T1_gff"] for s in SAMPLES]

# Import rules and functions
include: "rules/targets.smk"
include: "rules/snippy_pair.smk"
include: "rules/nucdiff.smk"
include: "rules/mlst.smk"
include: "rules/abricate.smk"
include: "rules/isescan.smk"
include: "rules/panaroo.smk"
include: "rules/phispy.smk"
include: "rules/sniffles.smk"

rule all:
    input:
        GhaisTargetFiles
