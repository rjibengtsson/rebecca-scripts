'''
Author: Rebecca

Description: This script will generate pairwise comparison of vcf files subsetted from multi vcf file

Usage:
'''

##################
# Import modules #
##################

import argparse
import subprocess
import pathlib

####################
# Take input files #
####################

parser = argparse.ArgumentParser(description='Generate pairwise comparison VCF file')
parser.add_argument('-vcf','--vcf', help='Merged VCF file', required=True)
parser.add_argument('-list','--list', help='List containing pairwise isolate accession', required=True)
args = parser.parse_args()

vcf_file = args.vcf
pw_list_file = args.list




def generate_pw_vcf_files(vcf_file, pw_list_file):

    patient = str(vcf_file).split("/")[0]
    path = pathlib.Path.cwd().joinpath(f'{patient}')
    pw_samples = open(pw_list_file, 'r').read().splitlines()

    for pair in pw_samples:
        pair = pair.split("\t")
        s1 = pair[0]
        s2 = pair[1]
        cmd = f'bcftools view -s {s1},{s2} --threads 4 -c 1  -Q 0.99  {vcf_file} > {path}/{s1}_vs_{s2}.vcf'
        subprocess.run([cmd], shell=True)

generate_pw_vcf_files(vcf_file, pw_list_file)