
#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


"""
get the non-synonymous/synonymous ratio for a set of genes (in nucletoide form)

requires biopython to be installed (conda install biopython)

"""


# Codon table with corresponding amino acids
codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': 'Stop', 'TAG': 'Stop', 'TGT': 'C', 'TGC': 'C', 'TGA': 'Stop', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def count_synonymous_non_synonymous_changes(codon):
    """
    calcs and returns the number of synonymous and non-synonymous changes for each codon 
    """
    synonymous_changes = 0
    non_synonymous_changes = 0
    
    original_aa = codon_table[codon]
    

    for i in range(3):
        for base in ['A', 'C', 'G', 'T']:
            mutated_codon = list(codon)
            mutated_codon[i] = base
            mutated_codon = ''.join(mutated_codon)
            
            if mutated_codon != codon:
                mutated_aa = codon_table[mutated_codon]
                if mutated_aa == original_aa:
                    synonymous_changes += 1
                else:
                    non_synonymous_changes += 1
    
    return synonymous_changes, non_synonymous_changes

# Calculate the ratio for each codon and store in a dictionary
codon_ratios = {}

for codon in codon_table.keys():
    synonymous_changes, non_synonymous_changes = count_synonymous_non_synonymous_changes(codon)
    
    if synonymous_changes == 0:
        ratio = 0 
    else:
        ratio = round(synonymous_changes / non_synonymous_changes, 2)

    
    codon_ratios[codon] = ratio

# Print the dictionary of codon ratios



def calculate_codon_usage(dna_sequences, codon_ratios):
    """
    get the synonymous/non-synonymous ratio for a sequence
    """

    codon_counts = {}
    total_codons = 0

    total_ratio_sum = 0

    for sequence in dna_sequences:
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i + 3].upper()
            if len(codon) == 3:  # Ignore any incomplete codons
                codon_counts[codon] = codon_counts.get(codon, 0) + 1
                total_codons += 1

                if codon in codon_ratios:
                    total_ratio_sum += codon_ratios[codon]


    return total_ratio_sum, total_codons


# arg parse to actually accept a file


def get_input():
	usage = 'python3 calc_codon_bias.py ...'
	parser = argparse.ArgumentParser(description='script to calculate crude non-synonymous/synonymous nucleotide change ratio for a set of genes. Requires biopython.', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file in fasta format (nucleotides)',  required=True)
	args = parser.parse_args()

	return args

args = get_input()

sample = args.infile.split('.')[0]

total_codons = 0
total_ratio_sum = 0

for dna_record in SeqIO.parse(args.infile, "fasta"):

    #Sample DNA sequence
    dna_sequence = dna_record.seq

    # Calculate the average ratio for the given sequence
    (ratio_sum, codon_count) = calculate_codon_usage([dna_sequence], codon_ratios)
    total_ratio_sum += ratio_sum
    total_codons += codon_count




average_ratio = total_ratio_sum / total_codons
average_ratio = 1/ average_ratio
print(total_codons)
print(f"Average Ratio: {average_ratio:.2f}")
    









