##################
# Import modules #
##################

import sys, os
import pandas as pd
import argparse
import re
import warnings
import numpy as np
warnings.simplefilter(action='ignore', category=FutureWarning)


####################
# Take input files #
####################

parser = argparse.ArgumentParser(description='Extract metadata from a csv file given a sample list')
parser.add_argument('meta', help='metadata file in txt format')
parser.add_argument('names', help='list of sample names')
parser.add_argument('out')
args = parser.parse_args()

meta_file = args.meta
sample_names = args.names
out = args.out

#############
# Functions #
#############

def extract_meta(meta_file, sample_names, out):
    sample_list = open(sample_names, "r").read().splitlines() # sample list
    df = pd.read_csv(meta_file, header=0, error_bad_lines=False, sep="\t") # metadata file


    col_header = list(df.columns)
    df2 = pd.DataFrame(columns = col_header)

    # convert column into list
    barcode_lst = list(df['SRA'])
    for sample in sample_list:
        # print(sample)
        if sample in barcode_lst:
            ## Get index from dataframe where the given value exists
            row_index = int(df[df['SRA'] == sample].index.values)
            ## select row from df
            data = df.iloc[[row_index]]
            df2 = df2.append(data)
        else:
            continue


    df2.to_csv(out, sep='\t', index=False)

########
# Main #
########

extract_meta(meta_file, sample_names, out)