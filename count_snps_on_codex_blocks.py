'''
Author: Rebecca

Description: This script will take a vcf file and count the number of SNPs identified on each codex block

Usage: python3 ~/scripts/python/count_snps_on_codex_blocks.py -vcf 0000.vcf -out snp_count.txt
'''

##################
# Import modules #
##################

import argparse
import subprocess

####################
# Take input files #
####################

parser = argparse.ArgumentParser(description='Count number of SNPs on codex blocks from a VCF file')
parser.add_argument('-vcf','--vcf', help="Input VCF file", required=True)
parser.add_argument('-out','--out', help="Output file", required=True)
args = parser.parse_args()

vcf_file = args.vcf
out_file = args.out


#############
# Functions #
#############


def removed_vcf_headers(vcf_file):

    in_file = open(vcf_file, 'r').readlines()

    with open('tmp.txt', 'w') as out_file:
        for line in in_file:
            if line.startswith("#"):
                pass
            else:
                out_file.write(line)

    return 'tmp.txt'


def create_codex_dict(vcf_file):

    in_file = open(removed_vcf_headers(vcf_file), 'r').readlines()
    codex_blocks_dict = {}

    for line in in_file:
        line = line.split('\t')[0]
        if line in codex_blocks_dict:
            counter += 1
        else:
            counter = 1
        codex_blocks_dict[line] = counter

    return codex_blocks_dict


def write_dict_to_table(vcf_file, out_file):

    dictionary = create_codex_dict(vcf_file)

    with open(out_file,'w') as out:
        for k, v in dictionary.items():
            line = k + "\t" + str(v) + "\n"
            out.write(line)

    tmp_file = removed_vcf_headers(vcf_file)
    cmd = f"rm {tmp_file}"
    subprocess.run([cmd], shell=True)

    return out_file


write_dict_to_table(vcf_file, out_file)