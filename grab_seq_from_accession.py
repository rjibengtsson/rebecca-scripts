'''
Description: This script will extract protein sequence given the name of the protein

Usage: python grab_seq_from_accession.py -acc NC_013364.1 -fasta GCF_000010765.1.fna -out GCF_000010765.1.chr.fna
'''


##################
# Import modules #
##################

import sys
from Bio import SeqIO
import re
import argparse

####################
# Take input files #
####################

parser = argparse.ArgumentParser(description='Extract nucleotide sequence from fna/faa file given the accession')
parser.add_argument('-acc', help='Accession list')
# parser.add_argument('-fasta', help='fna file name, must give absolute path of file')
# parser.add_argument('-out', help='name of output file')
args = parser.parse_args()

accession = args.acc
# fasta_file = args.fasta
# out_name = args.out

#############
# Functions #
#############


def find_seq(accession):
    acc_list = open(accession, 'r').read().splitlines()
    for acc in acc_list:
        refseq = acc.split("\t")[0]
        assembly_acc = acc.split("\t")[1]
        out_handle = open(refseq + ".chr.fna", 'w')
        print(assembly_acc)
        for record in SeqIO.parse(refseq + ".fna", 'fasta'):
            record_description = str(record.description)
            record_ID = record_description.split(" ")[0]
            if record_ID == assembly_acc:
                SeqIO.write(record, out_handle, 'fasta')

########
# Main #
########

find_seq(accession)
