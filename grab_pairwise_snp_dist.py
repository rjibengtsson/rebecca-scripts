'''
Description: This script will grab pairwise SNP distance from a sdiff_matrix.tsv file

Usage: python3 ~/scripts/python/grab_pairwise_snp_dist.py -s1 ERR2138465 -s2 ERR2138465 -tsv sdiff_matrix.tsv
'''

##################
# Import modules #
##################

import textwrap
import argparse
import pandas as pd
import os, sys

def parse_args():
    # set up command line arguments
    parser = argparse.ArgumentParser(
        prog = 'grab_pairwise_snp_dist.py',
        usage='%(prog)s [options] -tsv TSV',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        Script can be run two different ways:
        -------------------------------------
        Single mode: grab_pairwise_snp_dist.py -s1 S1 -s2 S2 -tsv TSV
        Batch mode: grab_pairwise_snp_dist.py -list LIST -tsv TSV
        '''),
        epilog = 'Get pairwise SNP distance from sdiff_matrix.tsv file.')
    parser.add_argument('-s1','--s1', help='sample 1 [optional]', type=str)
    parser.add_argument('-s2','--s2', help='sample 2 [optional]', type=str)
    parser.add_argument('-tsv','--tsv', help='sdiff_matrix.tsv file', required=True)
    parser.add_argument('-list','--list', help='list in text delimited format containing accessions of each pair [optional]')

    try:
        args = parser.parse_args()

        if args.list is None:
            if not (args.s1 and args.s2):
                sys.stderr.write(f"ERROR:Please provide {args.s1} and {args.s2}\n")
                sys.exit(1)

        if not os.path.exists(args.tsv):
            sys.stderr.write(f"ERROR:Please check input file path for {args.tsv}\n")
            sys.exit(1)

        if not args.list == None:
            if not os.path.exists(args.list):
                sys.stderr.write(f"ERROR:Please check input file path for {args.list}\n")
                sys.exit(1)

    except Exception as e:
        sys.stderr.write(str(e))
        sys.exit(1)

    return args



def read_sdiff_matrix(sdiff_file):

    sdiff_df = pd.read_csv(sdiff_file, sep="\t", header=0)
    sdiff_df2 = sdiff_df.set_index(sdiff_df.columns[0])

    return sdiff_df2


def get_pw_SNP_dist_single(sample1, sample2, sdiff_df):

    pairwise_snps = sdiff_df.loc[sample1, sample2]

    return print(f'{sample1} vs {sample2} \t {pairwise_snps}')


def get_pw_SNP_dist_loop(list_df, sdiff_df):


    with open('pairwise_snp_distance.txt', 'w') as outfile:
        for i in list_df.index:
            sample1 = list_df[0][i]
            sample2 = list_df[1][i]
            pairwise_snps = sdiff_df.loc[sample1, sample2]
            line = sample1 + "\t" + sample2 + "\t" + str(pairwise_snps) + "\n"
            outfile.write(line)

    return 'pairwise_snp_distance.txt'


def determine_single_or_loop(args):

    sdiff_file = args.tsv
    df = read_sdiff_matrix(sdiff_file)

    if args.list:
        list = pd.read_csv(args.list, sep="\t", header=None)
        get_pw_SNP_dist_loop(list, df)

    elif args.s1 and args.s2:
        sample1 = args.s1
        sample2 = args.s2
        get_pw_SNP_dist_single(sample1, sample2, df)


# Check if this file is run directly by python or is it being imported
if __name__ == "__main__":
    args = parse_args()
    determine_single_or_loop(args)

