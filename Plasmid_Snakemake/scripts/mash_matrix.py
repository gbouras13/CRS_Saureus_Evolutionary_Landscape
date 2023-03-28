#!/usr/bin/env python3

import pandas as pd


def parse_mash(input_mash_list, output_matrix):    

    mashes = []

    for mash_input in input_mash_list:
        tmp_mash = read_mash(mash_input)
        mashes.append(tmp_mash)

    # make into combined dataframe
    total_mash_df = pd.concat(mashes,  ignore_index=True)

    # get only sample name
    total_mash_df['query'] = total_mash_df['query'].str.replace('.fasta','')
    #https://stackoverflow.com/questions/15012228/splitting-on-last-delimiter-in-python-string
    total_mash_df[['query', 'query_sample']] = total_mash_df['query'].str.rsplit('/', 1, expand=True)
    #total_mash_df['query'] = total_mash_df['query'].str.rsplit('/', 1)[1]

    total_mash_df['reference'] = total_mash_df['reference'].str.replace('.fasta','')
    total_mash_df[['reference', 'reference_sample']] = total_mash_df['reference'].str.rsplit('/', 1, expand=True)


    total_mash_pivoted = total_mash_df.pivot_table(index='query_sample', columns='reference_sample', values='distance')

    
    total_mash_pivoted.to_csv(output_matrix, sep=",", index=False)


def read_mash(mash_input):
    colnames=['query', 'reference', 'distance', 'p_val', 'hashes'] 
    mash_df = pd.read_csv(mash_input, delimiter= '\t',  header=None, names=colnames)
    return(mash_df)
                
parse_mash(snakemake.input.mashes,snakemake.output.matrix)



        
