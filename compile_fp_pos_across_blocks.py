'''
Author: Rebecca

Description: This script will take the concatenated false positive vcf file and count how many strains
has false positive SNPs at a particular site.
'''

import sys

vcf_file = open(sys.argv[1], "r")
tmp_file = open('tmp.txt' , 'w')

# reformat file to remove lines containing headers


