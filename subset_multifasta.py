'''
This script will subset a multifasta according to gene names
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

parser = argparse.ArgumentParser(description='Extract gene sequences from multifasta file given the name of the gene')
parser.add_argument('name', help='list containing gene names')
parser.add_argument('fasta', help='multifasta file you want to subset')
args = parser.parse_args()

gene_name = args.name
fasta_file = args.fasta

#############
# Functions #
#############

def subset_fasta(gene_name, fasta_file):

    # convert gene names into a list
    gene_lst = list(open(gene_name).read().splitlines())

    for gene in gene_lst:

        # create a new multifasta file for each gene
        out_handle = open(gene + '.fasta', 'w')

        # loop through the multifasta file
        for record in SeqIO.parse(fasta_file, 'fasta'):
            record_name = str(record.id).split(":")[0]
            if record_name == gene:
                SeqIO.write(record, out_handle, 'fasta')
            else:
                continue

########
# Main #
########

subset_fasta(gene_name, fasta_file)