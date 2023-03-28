#!/usr/bin/env python3

import pandas as pd
import numpy as np 

def parse_jaccard(input_rtab, output_matrix):    

    rtab_matrix = pd.read_csv(input_rtab, delimiter= '\t', index_col=False, header=0)
    # remove first column (genes)
    rtab_matrix = rtab_matrix.iloc[: , 1:]

    # instantiate the matrix
    jac_mat = np.zeros((len(rtab_matrix.columns), len(rtab_matrix.columns)))

    for i in range(len(rtab_matrix.columns)):
        for j in range(len(rtab_matrix.columns)):
            df_tmp = rtab_matrix.iloc[:, [i, j]]
            # drop 0s https://stackoverflow.com/questions/22649693/drop-rows-with-all-zeros-in-pandas-data-frame
            df_tmp = df_tmp.loc[~(df_tmp==0).all(axis=1)]
            df_tmp = df_tmp.set_axis(['A', 'B'], axis=1, inplace=False)
            # if both share gene, then sum will be 2
            df_tmp['sum'] = df_tmp['A'] + df_tmp['B']
            counts = (df_tmp['sum'] == 2).sum()
            jac_mat[i][j] = counts/len(df_tmp.index)

    jac_mat = pd.DataFrame(jac_mat, index=rtab_matrix.columns, columns=rtab_matrix.columns)
   
    jac_mat.to_csv(output_matrix, sep=",", index=False)


                
parse_jaccard(snakemake.input.rtab, snakemake.output.matrix)



        
