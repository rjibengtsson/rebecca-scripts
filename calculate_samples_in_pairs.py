'''
Description: This script will calculate number of samples contributing to pairs from the Ludden study

Usage: python3 ~/scripts/python/calculate_samples_in_pairs.py -input Gpx_vs_Ludden_09_11_2022_diff.csv -output accessionID1_occurence.txt
'''

##################
# Import modules #
##################

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Calculate how many times individual samples contribute to pairs')
parser.add_argument('-input','--input', help='Input file containing pairwise accession in column 1 and 2', required=True)
parser.add_argument('-output','--output', help='Output file', required=True)

args = parser.parse_args()

input_file = args.input
output_file = args.output


def remove_recombination_pairs(input_file):

    df = pd.read_csv(input_file, header=0)
    # print(len(df.index))
    df2 = df[df['SNPs after removal of recombination '] < 100]
    # print(len(df2.index))

    return df2


def count_isolate_occurence_in_pairs(input_file):

    df = remove_recombination_pairs(input_file)

    acc_1_list = list(df['AccessionID2'])

    acc_dict = {}

    for acc in acc_1_list:
        if acc in acc_dict:
            acc_dict[acc] +=1
        else:
            acc_dict[acc] = 1

    acc_dict_srt = sorted(acc_dict.items(), key=lambda x:x[1], reverse=True)
    converted_acc_dict_srt = dict(acc_dict_srt)

    return converted_acc_dict_srt


def write_dict_to_table(input_file, output_file):

    acc_dict = count_isolate_occurence_in_pairs(input_file)

    with open(output_file, 'w') as out:
        for k, v in acc_dict.items():
            line = k + "\t" + str(v) + "\n"
            out.write(line)


write_dict_to_table(input_file, output_file)