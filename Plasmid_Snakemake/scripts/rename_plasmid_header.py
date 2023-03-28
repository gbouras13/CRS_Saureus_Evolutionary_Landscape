#!/usr/bin/env python3

from Bio import SeqIO

# get chromosome if length > 2.5 million

def rename_header( input_fasta, output_fasta, plasmid):
    # read in the fasta
    with open(output_fasta, 'w') as fa:
        for dna_record in SeqIO.parse(input_fasta, "fasta"):
            dna_record.id = plasmid
            dna_record.description = ""
            SeqIO.write(dna_record, fa, 'fasta')

            
rename_header(snakemake.input[0],snakemake.output[0], snakemake.wildcards.plasmid)



        
