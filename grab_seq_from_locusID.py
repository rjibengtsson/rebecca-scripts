### This script will extract protein sequence given the name of the protein ###


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

parser = argparse.ArgumentParser(description='Extract protein sequence from faa file given the name of the protein')
parser.add_argument('-name', help='protein name')
parser.add_argument('-file', help='fna file name, must give absolute path of file')
parser.add_argument('-out', help='name of output file')
args = parser.parse_args()

protein_name = args.name
file_name = args.file
out_name = args.out

#############
# Functions #
#############


def find_protein(protein_name, file_name, out_name):
    protein = open(protein_name).read().splitlines()
    out_handle = open(out_name, 'a')
    for i in protein:
        match = '[locus_tag=' + i + ']'
        for record in SeqIO.parse(file_name, 'fasta'):
            record_description = str(record.description)
            record_ID = record_description.split(" ")[1]
            if record_ID == match:
                SeqIO.write(record, out_handle, 'fasta')

########
# Main #
########

find_protein(protein_name, file_name, out_name)