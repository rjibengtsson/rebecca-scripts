'''
Author: Rebecca

Description: This script will take a list of names and extract relevant info from a metadata file

Usage: python3
'''


##################
# Import modules #
##################

import sys, os
import textwrap
import pandas as pd
import argparse
import re
import warnings
import numpy as np
warnings.simplefilter(action='ignore', category=FutureWarning)


####################
# Take input files #
####################

def parse_args():
    # set up command line arguments
    parser = argparse.ArgumentParser(
        prog = extract_metadata.py,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
            Subset data from a metadata file
            '''))
    parser.add_argument('-meta', '--meta', help='metadata file', required=True)
    parser.add_argument('-list', '--list', help='list of samples to subset', required=True)
    parser.add_argument('-col', '--col', help='column name to match against your list', required=True)
    parser.add_argument('-fmt', '--fmt', help='metadata file format, txt or csv [Default: txt]', default="txt", type=str)
    parser.add_argument('-out', '--out', help='output file', required=True)

    try:
        args = parser.parse_args()

        if not os.path.exists(args.meta):
            sys.stderr.write(f"ERROR:Please check input file path for {args.meta}\n")
            sys.exit(1)

        if not os.path.exists(args.list):
            sys.stderr.write(f"ERROR:Please check input file path for {args.list}\n")
            sys.exit(1)

    except Exception as e:
        sys.stderr.write(str(e))
        sys.exit(1)

    return args

#############
# Functions #
#############

def read_metadata_to_df(args):

    # Determine if metadata file provided is csv or txt format
    if args.fmt == "csv":
        df = pd.read_csv(args.meta, header=0, error_bad_lines=False)
    else:
        df = pd.read_csv(args.meta, header=0, error_bad_lines=False, sep="\t")

    return df


def extract_meta(args):

    df = read_metadata_to_df(args)
    col_header = list(df.columns)
    df2 = pd.DataFrame(columns=col_header)

    sample_list = open(sample_names, "r").read().splitlines() # sample list

    # convert column into list

    barcode_lst = list(df['Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Status)'])

    sra_acc_lst = []

    for i in barcode_lst:
        acc = i.split(";")[0]
        sra_acc_lst.append(acc)

    df['Accession No'] = sra_acc_lst

    sample_id_lst = list(df['Sample ID'])

    for sample in sample_list:
        # print(sample)
        if sample in sample_id_lst:
            ## Get index from dataframe where the given value exists
            row_index = int(df[df['Sample ID'] == sample].index.values)
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