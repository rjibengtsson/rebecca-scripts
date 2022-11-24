'''
Author: Rebecca

Description: This script will perform blastn and find the coordinates of key genes
from phase I input reference strains
'''

##################
# Import modules #
##################

import argparse
import subprocess
from Bio import SeqIO


####################
# Take input files #
####################


parser = argparse.ArgumentParser(description='find key gene coordinates in input reference strains')
parser.add_argument('-query','--query', help='key gene fasta file', required=True)
parser.add_argument('-subject','--subject', help='phase I strains fasta file', required=True)
parser.add_argument('-gff','--gff', help='phase I strains gff file', required=True)
args = parser.parse_args()

query_fasta = args.query
subject_fasta = args.subject
gff_file = args.gff


# set key gene name
key_gene_name = str(query_fasta).split('.')[0]


#############
# Functions #
#############


# perform makeblastbd on input reference strains
def make_blastndb(subject_fasta):
    cmd = f'makeblastdb -in {subject_fasta} -input_type fasta -dbtype nucl -parse_seqids'
    blastdb = subprocess.run([cmd], shell=True)
    return blastdb


# blast key gene sequence against reference strain database
def blastn_keyGene(subject_fasta, query_fasta):
    cmd = f'blastn -db {subject_fasta} -query {query_fasta} ' \
          f'-outfmt 6 -out {key_gene_name}_vs_phaseI_ref.raw.blastn'
    subprocess.run([cmd], shell=True)
    return f'{key_gene_name}_vs_phaseI_ref.raw.blastn'


# calculate key gene length
def calc_keyGene_length(query_fasta):
    record = SeqIO.read(query_fasta, "fasta")
    record_length = len(record.seq)
    return int(record_length)


# filter blastn output, grabbing only hits with >99% sequence coverage
def filter_blastn_results():
    query_length = calc_keyGene_length(query_fasta)
    raw_blast_file = blastn_keyGene(subject_fasta, query_fasta)

    with open(f'{key_gene_name}_vs_phaseI_ref.filt.blastn', 'w') as tmp_file:
        with open(raw_blast_file, 'r') as file:
            for line in file:
                index = line.split()
                blastn_len = int(index[3])
                seq_coverage = round(blastn_len / query_length, 2)
                if seq_coverage >= 0.99:
                    tmp_file.write(line)
                else:
                    pass
    return f'{key_gene_name}_vs_phaseI_ref.filt.blastn'


# turn blastn file into bed file
def make_bed_file():
    blastn_file = filter_blastn_results()
    bed_file = blastn_file + '.bed'

    with open(bed_file, 'w') as out_file:
        with open(blastn_file, 'r') as file:
            for line in file:
                line = line.split()
                new_line = line[1] + '\t' + line[8] + '\t' + line[9] + '\t' + line[0] + '\n'
                out_file.write(new_line)
    return bed_file


# sort the bed file so that it's compatible to run with bedtools intersect
def sort_bed_file():
    bed_file = make_bed_file()
    sorted_bed_file = str(bed_file).split('.')[0] + '.srt.bed'

    with open('tmp.file', 'w') as tmp_file:
        with open(bed_file, 'r') as infile:
            in_file = infile.readlines()
            for line in in_file:
                index = line.split()
                if int(index[1]) > int(index[2]):
                    new_line = index[0] + '\t' + index[2] + '\t' + index[1] + '\t' + index[3] + '\n'
                    tmp_file.write(new_line)
                else:
                    tmp_file.write(line)

    cmd = 'sort -k1,1 -k2,2n tmp.file > ' + sorted_bed_file
    subprocess.run([cmd], shell=True)

    cmd2 = 'rm tmp.file'
    subprocess.run([cmd2], shell=True)

    return sorted_bed_file


# format original gff file and display only pgap annotated genes
def format_gff(gff_file):
    pgap_gff_file = str(gff_file).split(".")[0] + '.pgap.gff'

    with open(pgap_gff_file, 'w') as out_file:
        with open(gff_file, 'r') as file:
            for line in file:
                if not line.startswith("#"):
                    index = line.split()
                    if index[2] == 'gene':
                        new_line = index[0] + "\t" + index[3] + "\t" + index[4] + "\t" + index[8] + "\n"
                        out_file.write(new_line)
                else:
                    continue
    return pgap_gff_file


# perform bedtools intersect
def bedtools_intersect():

    gff = format_gff(gff_file)
    bedfile = sort_bed_file()
    bed_intersect_file = f'{key_gene_name}_vs_phaseI_ref.bedintersect'

    cmd = f'bedtools intersect -a {gff} -b {bedfile} -wa > {bed_intersect_file}'
    subprocess.run([cmd], shell=True)
    return


########
# Main #
########


make_blastndb(subject_fasta)
bedtools_intersect()