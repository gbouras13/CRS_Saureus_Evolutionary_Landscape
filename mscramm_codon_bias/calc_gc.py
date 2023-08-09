
#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


"""
get the overall GC ratio 

"""



def get_input():
	usage = 'python3 calc_gc.py ...'
	parser = argparse.ArgumentParser(description='script to calculate crude non-synonymous/synonymous nucleotide change ratio for a set of genes. Requires biopython.', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file in fasta format (nucleotides)',  required=True)
	args = parser.parse_args()

	return args

args = get_input()


sequences = SeqIO.parse(args.infile, "fasta")


total_gc_count = 0
total_base_count = 0
    
for record in sequences:
    sequence = str(record.seq)
    gc_count = sequence.count('G') + sequence.count('C')
    base_count = len(sequence)
    
    total_gc_count += gc_count
    total_base_count += base_count

overall_gc_percentage = (total_gc_count / total_base_count) * 100
print(f"Overall GC Percentage: {overall_gc_percentage:.2f}%")
    









