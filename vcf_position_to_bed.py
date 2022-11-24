''''
Description: This script will take a VCF file and convert the positions into a bed file

Usage:
'''


import argparse

####################
# Take input files #
####################

parser = argparse.ArgumentParser(description='Convert variant positions to bed positions')
parser.add_argument('-vcf','--vcf', help="Input VCF file", required=True)
parser.add_argument('-out','--out', help="Output bed file", required=True)
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

def convert_to_bed(vcf_file, out_file):

    in_file = open(removed_vcf_headers(vcf_file), 'r').readlines()

    with open(out_file, 'w') as out:
        for line in in_file:
            line = line.split('\t')
            block_name = line[0]
            start = int(line[1]) - 1
            end = line[1]
            new_line = block_name + "\t" + str(start) + "\t" + str(end) + "\n"
            out.write(new_line)

    return out_file

convert_to_bed(vcf_file, out_file)