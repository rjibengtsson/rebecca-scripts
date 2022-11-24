'''
Author: Rebecca

Description: This script will perform blastn key gene against codex and grab only hits with >99% sequence coverage
'''

import argparse
import subprocess
from Bio import SeqIO


parser = argparse.ArgumentParser(description='blastn key gene against codex')
parser.add_argument('-query','--query', help='key gene fasta file', required=True)
parser.add_argument('-subject','--subject', help='codex fasta file', required=True)
args = parser.parse_args()

query_fasta = args.query
subject_fasta = args.subject

key_gene_name = str(query_fasta).split('.')[0]


def make_blastndb(subject_fasta):
    cmd = f'makeblastdb -in {subject_fasta} -input_type fasta -dbtype nucl -parse_seqids'

    blastdb = subprocess.run([cmd], shell=True)

    return blastdb


def blast_keygene(subject_fasta, query_fasta):

    cmd = f'blastn -db {subject_fasta} -query {query_fasta} ' \
          f'-outfmt 6 -out {key_gene_name}.raw.blastn'

    subprocess.run([cmd], shell=True)

    return f'{key_gene_name}.raw.blastn'


def calc_keyGene_length(query_fasta):
    record = SeqIO.read(query_fasta, "fasta")
    record_length = len(record.seq)
    return int(record_length)


def filter_blastn_results():
    query_length = calc_keyGene_length(query_fasta)

    # filter blastn output
    raw_blast_file = blast_keygene(subject_fasta, query_fasta)

    with open(f'{key_gene_name}.filt.blastn', 'w') as tmp_file:
        with open(raw_blast_file, 'r') as file:
            for line in file:
                index = line.split()
                blastn_len = int(index[3])
                seq_coverage = round(blastn_len / query_length, 2)
                if seq_coverage >= 0.99:
                    tmp_file.write(line)
                else:
                    pass

    return f'{key_gene_name}.filt.blastn'




make_blastndb(subject_fasta)
filter_blastn_results()