'''
Author: Rebecca

Description: This script will take a list of names and extract relevant info from a metadata file

Usage: python ~/scripts/python/extract_metadata.py -meta enterobase_ecoli_clermontype_meta.txt -list low_res_samples.txt -col "SRA" -fmt
txt -out low_res_samples_clermontType.txt
'''


##################
# Import modules #
##################

import sys, os
import textwrap
import pandas as pd
import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


####################
# Take input files #
####################

def parse_args():
    # set up command line arguments
    parser = argparse.ArgumentParser(
        prog = 'extract_metadata.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
            Subset data from a metadata file
            '''))
    parser.add_argument('-meta', '--meta', help='metadata file', required=True)
    parser.add_argument('-list', '--list', help='list of samples to subset', required=True)
    parser.add_argument('-col', '--col', help='column name to match against your list, must be provided within double qoutes', required=True)
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

    sample_list = open(args.list, "r").read().splitlines() # sample list

    sample_id_lst = list(df[args.col])

    for sample in sample_list:
        if sample in sample_id_lst:
            ## Get index from dataframe where the given value exists
            row_index = int(df[df[args.col] == sample].index.values)
            ## select row from df
            data = df.iloc[[row_index]]
            df2 = df2.append(data)
        else:
            empty_dict = {}
            for element in col_header:
                if element == args.col:
                    empty_dict[element] = [sample]
                else:
                    empty_dict[element] = ["NA"]
            empty_df = pd.DataFrame(empty_dict)
            df2 = pd.concat([df2, empty_df])


    df2.to_csv(args.out, sep='\t', index=False)


########
# Main #
########

if __name__ == "__main__":
    args = parse_args()

    # Call functions
    extract_meta(args)