'''
This script will grab gene length information from a gff file
'''

import sys
import pandas as pd

# read reformated gff file into dataframe
df = pd.read_csv(sys.argv[1], header=None, sep="\t")

gene_len_lst = []

# iterate through dataframe
for index, row in df.iterrows():
    gene_len = row[4] - row[3]
    gene_len_lst.append(gene_len)

# append list as column in dataframe
df['gene_len'] = gene_len_lst

# write new dataframe into csv
df.to_csv(sys.argv[2], index=False)
