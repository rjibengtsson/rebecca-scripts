'''
Description: This script will grab pairwise SNP distance from a sdiff_matrix.tsv file

Usage: python3 ~/software/scripts/grab_pairwise_snp_dist_from_list.py -list pw_acc_list.txt -tsv sdiff_matrix.tsv -out pw_snp_dist.txt
'''

##################
# Import modules #
##################

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Get pairwise SNP distance for two samples')
parser.add_argument('-list','--list', help='List containing pairwise sample accessions', required=True)
parser.add_argument('-tsv','--tsv', help='sdiff_matrix.tsv file', required=True)
parser.add_argument('-out','--out', help='Output file', required=True)
args = parser.parse_args()

sample_lst = args.list
sdiff_matrix = args.tsv
out_file = args.out

def get_pw_SNP_dist(sample1, sample2, sdiff_matrix):

    df = pd.read_csv(sdiff_matrix, sep="\t", header=0)
    df2 = df.set_index(df.columns[0])

    pairwise_snps = df2.loc[sample1, sample2]
    line = sample1 + "\t" + sample2 + "\t" + str(pairwise_snps) + "\n"
    return line

def read_pairwise_list(sample_lst, out_file, sdiff_matrix):

    pw_lst = open(sample_lst, 'r').read().splitlines()

    with open(out_file, 'w') as out:
        for pw in pw_lst:
            sample1 = pw.split()[0]
            sample2 = pw.split()[1]
            line = get_pw_SNP_dist(sample1, sample2, sdiff_matrix)
            out.write(line)

    return out_file


read_pairwise_list(sample_lst, out_file, sdiff_matrix)