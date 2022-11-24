### This script will take a tab delimited file and change the name accordingly ###

##################
# Import modules #
##################

import os, sys
import argparse


parser = argparse.ArgumentParser(description='Change file name')
parser.add_argument('-f', help='tab delimited file containing sample names')
# parser.add_argument('-opath', help='path of output files, insert with double quote')
args = parser.parse_args()

in_file = args.f
# out_path = args.opath

#############
# Functions #
#############

# add dates to end of fasta header
def change_name(in_file):
    name_file = open(in_file, 'r').read().splitlines()
    for name in name_file:
        acc = name.split("\t")[1]
        refseq = name.split("\t")[0]
        cmd = "mv " + refseq + ".fastq.gz  " + acc + ".fastq.gz "
        print(cmd)
        # os.system(cmd)

########
# Main #
########

change_name(in_file)