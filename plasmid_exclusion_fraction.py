'''
Description: This script will grab the total size of noplasmid.fasta file and the Screened_Unitig_Size,
calculate what percentage of the genome is excluded beacuse of plasmid exclusion

Usage: python3 plasmid_exclusion_fraction.py -list samples.txt -out plasmid_exclusion_fraction.txt
'''

import argparse
import pandas as pd
from Bio import SeqIO


parser = argparse.ArgumentParser(description='Calculate percentage of genome lost due to plasmid exclusion')
parser.add_argument('-list','--list', help='List containing sample names', required=True)
parser.add_argument('-out','--out', help='Out file containing final output', required=True)
args = parser.parse_args()

sample_lst = args.list
out_file = args.out

def grab_unitig_size_from_QC(sample_lst):

    samples = open(sample_lst, 'r').read().splitlines()

    size_dict = {}

    for sample in samples:
        file = f"./{sample}/QC.log"
        df = pd.read_csv(file, header=None, sep="\t")
        Assembly_Unitig_Size = df.loc[9, 1]
        noPlasmid_Unitig_Size = df.loc[13, 1]
        size_dict[sample] = [Assembly_Unitig_Size, noPlasmid_Unitig_Size]

    return size_dict


def calculate_plasmid_exclusion_size():

    size_dict = grab_unitig_size_from_QC(sample_lst)

    for key, value in size_dict.items():
        fraction = round(1 - (float(value[1]) / float(value[0])), 2)
        size_dict[key].append(fraction)

    return size_dict


def write_dict_to_table(out_file):

    size_dict = calculate_plasmid_exclusion_size()

    with open(out_file, 'w') as out:
        title = "Sample" + "\t" + "Unitigs assembly size" + "\t" + "Plasmid excluded unitigs assembly size" + "\t" + "Fraction excluded" + "\n"
        out.write(title)
        for key, value in size_dict.items():
            line = key + "\t" + str(value[0]) + "\t" + str(value[1]) + "\t" + str(value[2]) + "\n"
            out.write(line)

write_dict_to_table(out_file)