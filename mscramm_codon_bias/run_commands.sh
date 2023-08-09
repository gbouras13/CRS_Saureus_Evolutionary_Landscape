#!/bin/bash

# assumes biopython is installed

# to split the pan genome into MSCRAMM and non-MSCRAMM genes

python extract_mscramm_and_not_mscramm.py -i pan_genome_reference.fa

# to calculate the gc contents 

python calc_gc.py -i mscramm.fasta

python calc_gc.py -i non_mscramm.fasta

# to calculate the codon usage bias (ratio of non-synonymous to synonymous nucleotide changes for all codons)

python calc_codon_bias.py  -i mscramm.fasta

python calc_codon_bias.py -i non_mscramm.fasta