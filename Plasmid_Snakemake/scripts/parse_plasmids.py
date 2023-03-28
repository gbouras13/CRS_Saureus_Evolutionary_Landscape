#!/usr/bin/env python3

from plsdbapi import query
import pandas as pd
import os
from Bio import SeqIO


def parse_plasmids(input_fasta, tmp_dir,  output_dir, output_touch, sample):    

    for dna_record in SeqIO.parse(input_fasta, "fasta"):
        contig_file = os.path.join(tmp_dir, sample + "_" + dna_record.id + ".fasta")
        with open(contig_file, 'w') as fa:
            SeqIO.write(dna_record, fa, 'fasta')
        #print(contig_file)
        # error handling if the query fails
        try:
            df = query.query_plasmid_sequence('mash_screen', ifile=contig_file, mash_max_v=1, mash_min_i=0)
            df.to_csv(os.path.join(output_dir, sample + "_" + dna_record.id + ".csv"), index=False) 
        except:
            touch(os.path.join(output_dir, sample + "_" + dna_record.id + ".csv"))
            continue
    
    # write regardless if the file is empty for the snakemake rule
    touch(output_touch)
        

# function to touch create a file 
# https://stackoverflow.com/questions/12654772/create-empty-file-using-python
def touch(path):
    with open(path, 'a'):
        os.utime(path, None)
                
parse_plasmids(snakemake.input[0],snakemake.params[0], snakemake.params[1], snakemake.output[0], snakemake.wildcards.sample)



        
