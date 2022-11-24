'''
This script will sort bed files
'''

import sys, os
import subprocess


# reverse minus coordinates
def reverse_minus_coord(in_file):
    with open('tmp.bed', 'w') as tmp_file:
        with open(in_file, 'r') as infile:
            in_file = infile.readlines()
            for line in in_file:
                line = line.split()
                if int(line[1]) > int(line[2]):
                    new_line = line[0] + '\t' + line[2] + '\t' + line[1] + '\t' + line[3] + '\n'
                    tmp_file.write(new_line)
                else:
                    new_line = line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\n'
                    tmp_file.write(new_line)
    return 'tmp.bed'

# sort bed file
def sort_bed_file(infile, out_file):
    tmp_file = reverse_minus_coord(infile)
    cmd = 'sort -k1,1 -k2,2n ' + tmp_file + ' > ' + out_file
    result = subprocess.run([cmd], shell=True)
    return result


infile = sys.argv[1]
out_file = sys.argv[2]

sort_bed_file(infile, out_file)