'''
Author: Rebecca

Description: This script will take Ben's CSV file and grab all the false positive positions from the codex blocks,
then use a R script to visualise the distribution of the positions on each block.
'''

##################
# Import modules #
##################

import argparse
import pandas as pd
import subprocess

####################
# Take input files #
####################

parser = argparse.ArgumentParser(description='Display false positive positions on codex blocks')
parser.add_argument('-csv','--csv', help="Ben's csv file", required=True)
args = parser.parse_args()

csv_file = args.csv


#############
# Functions #
#############


def grab_FP_rows(csv_file):

    in_file = open(csv_file, 'r').readlines()
    out_file_name = str(csv_file).split(".")[0]

    with open(f'{out_file_name}_FPs.csv', 'w') as out_file:
        for line in in_file:
            if "#" in line or "FP" in line:
                out_file.write(line)
            else:
                pass

    return f'{out_file_name}_FPs.csv'


def create_FPs_dicts():

    in_file = open(grab_FP_rows(csv_file), 'r').readlines()
    codex_blocks_dict = {}

    for line in in_file[1:]:
        index = line.strip('\n').split(",")
        fp_pos = index[1]
        strain_fp_num = list(index[5:]).count("FP")
        fp_pos_dict = {fp_pos : strain_fp_num}
        if index[0] in codex_blocks_dict:
            codex_blocks_dict[index[0]][fp_pos] = strain_fp_num
        else:
            codex_blocks_dict[index[0]] = fp_pos_dict

    return codex_blocks_dict


def R_plot(in_file, codex_block_name):

    R_cmd = f"Rscript plot_FP_snp_distribution.R {in_file} {codex_block_name}.png "
    return subprocess.run([R_cmd], shell=True)



def format_dataframe_from_dict():

    dict = create_FPs_dicts()

    for key, value in dict.items():
        block_name = str(key)
        block_length = block_name.split("-")[3].replace("L", "")
        with open(f'{block_name}_df.txt', 'w') as out_file:
            for k in value.keys():
                position = k
                number_of_strains = value[k]
                line = block_name + "\t" + str(block_length) + "\t" + str(position) + "\t" + str(number_of_strains) + "\n"
                out_file.write(line)


def generate_R_plots():

    dict = create_FPs_dicts()

    for key, value in dict.items():
        codex_block_name = str(key)
        R_plot(f"{codex_block_name}_df.txt", codex_block_name)


########
# Main #
########


format_dataframe_from_dict()
generate_R_plots()