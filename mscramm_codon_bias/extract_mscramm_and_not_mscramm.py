#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def get_input():
	usage = 'python3 extract_mscramm_and_not_mscramm.py ...'
	parser = argparse.ArgumentParser(description='script to extract mscramm and non-mscramm genes. Requires biopython.', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file in fasta format (nucleotides)',  required=True)
	args = parser.parse_args()

	return args

args = get_input()



mscramm_genes = ["isdB",
"clfB",
"mscL",
"isdF",
"isdD",
"isdE",
"ebhB~~~ebh~~~ebhA",
"efb",
"isdA",
"sdrC",
"clfA",
"isdC",
"eap",
"fnbA",
"sdrD~~~clfA",
"sdrE",
"fnbB",
"cna",
"sasG"]

outfile = "mscramm.fasta"

with open(outfile, 'w') as out_fa:
    for dna_record in SeqIO.parse(args.infile, "fasta"):
        
        if dna_record.id in mscramm_genes:
            SeqIO.write(dna_record, out_fa, 'fasta')

outfile = "non_mscramm.fasta"

with open(outfile, 'w') as out_fa:
    for dna_record in SeqIO.parse(args.infile, "fasta"):
        
        if dna_record.id not in mscramm_genes:
            SeqIO.write(dna_record, out_fa, 'fasta')