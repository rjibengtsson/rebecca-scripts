'''
This script will subset isolates from a molten pairwise dist file
'''

import sys
import argparse
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


parser = argparse.ArgumentParser(description='Extract information from a csv file given a sample list')
parser.add_argument('meta', help='metadata file in csv format')
parser.add_argument('names', help='list of sample names')
parser.add_argument('out')
args = parser.parse_args()

meta_file = args.meta
sample_names = args.names
out = args.out

#############
# Functions #
#############

def data_extract(sample_names, meta_file, out):
    isolates_lst = list(open(sample_names, 'r').read().splitlines())

    df = pd.read_csv(meta_file, header=0)

    col_header = list(df.columns)
    df2 = pd.DataFrame(columns = col_header)

    for i in range(len(df)):
        if (df.loc[i, "col"] in isolates_lst) & (df.loc[i, "row"] in isolates_lst):
            data = df.iloc[i]
            df2 = df2.append(data)

    df2.to_csv(out, sep='\t', index=False)

data_extract(sample_names, meta_file, out)