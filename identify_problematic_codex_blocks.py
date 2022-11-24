'''
Description: This script will count how many times discordant pairwise SNPs occur and which codex block it occurs on,
then tallies up the how many times the same blocks are causing the discrepancy.

Usage: python identify_problematic_codex_blocks.py -list patient_vcf_file_list.txt -out problematic_blocks
'''

##################
# Import modules #
##################

import argparse
import pandas as pd
import pathlib

parser = argparse.ArgumentParser(description='Identify problematic codex blocks')
parser.add_argument('-list','--list', help='List containing the name of pairwise comparison vcf files', required=True)
parser.add_argument('-output','--output', help='Output file prefix', required=True)

args = parser.parse_args()

list_file = args.list
output_file = args.output



def create_codex_dict(vcf_file_path, patient_id, vcf_file_name):

    # create dictionary
    in_file = open(vcf_file_path, 'r').readlines()
    codex_blocks_dict = {}

    for line in in_file:
        if not line.startswith("#"):
            line = line.split('\t')[0]
            if line in codex_blocks_dict:
                counter += 1
            else:
                counter = 1
            codex_blocks_dict[line] = [counter]

    # append patient id and vcf file name to dictionary

    for key, value in codex_blocks_dict.items():
        codex_blocks_dict[key].append(patient_id)
        codex_blocks_dict[key].append(vcf_file_name)

    return codex_blocks_dict


def write_dict_to_dataframe(list_file):

    vcf_file_list = open(list_file, 'r').read().splitlines()
    df = pd.DataFrame(columns=['Codex block name', 'Num of discrepant SNPs', 'Patient ID', 'vcf file name'])

    for i in vcf_file_list:
        i = i.split('\t')
        vcf_file_name =i[1]
        patient_id = i[0]
        path = pathlib.Path.cwd().joinpath(f'{i[0]}/{vcf_file_name}')
        block_dict = create_codex_dict(path,patient_id,vcf_file_name)
        for k, v in block_dict.items():
            append_list = []
            append_list.append(k)
            for i in v:
                append_list.append(i)
            df.loc[len(df)] = append_list

    return df

def write_df_to_table(output_file):

    df = write_dict_to_dataframe(list_file)
    df.to_csv(output_file + '.txt', index=False, sep='\t')

    return output_file + '.txt'

def count_problem_blocks(output_file):

    file = write_df_to_table(output_file)
    df = pd.read_csv(file, header=0, sep="\t")
    patient_dict = {}

    for i in range(len(df)):
        patient_id = df.loc[i, "Patient ID"]
        block_name = df.loc[i, "Codex block name"]
        num_snps = int(df.loc[i, "Num of discrepant SNPs"])
        # If patient already exist in dictionary, iterate through block_dict
        if patient_id in patient_dict:
            # If block name (key) already exist, add to block_counter and sum number of snps
            if block_name in block_dict.keys():
                current_counter , current_num_snps = block_dict[block_name]
                new_counter = current_counter + 1
                new_num_snps = current_num_snps + num_snps
                block_dict[block_name] = [new_counter, new_num_snps]
            # Otherwise create new key for new block
            else:
                new_block_dict = {}
                block_counter = 1
                new_block_dict[block_name] = [block_counter, num_snps]
                patient_dict[patient_id].update(new_block_dict)
        # Otherwise
        else:
            # If encounters a new patient ID, create new block dictionary
            block_dict = {}
            block_counter = 1
            # In the block_dict set block_name to be key
            block_dict[block_name] = [block_counter, num_snps]
            patient_dict[patient_id] = block_dict

    return patient_dict

def write_patient_dict_to_file(output_file):

    patient_dict = count_problem_blocks(output_file)
    df = pd.DataFrame(columns=['Patient ID', 'Codex block name', 'Total recurrence','Total num of discrepant SNPs',])

    for key, value in patient_dict.items():
        patient_id = key
        for k, v in value.items():
            append_list = []
            append_list.append(patient_id)
            codex_block = k
            recurrence, total_num_snps = v
            append_list.append(codex_block)
            append_list.append(recurrence)
            append_list.append(total_num_snps)
            df.loc[len(df)] = append_list

    df.to_csv(output_file + '_counter.txt', index=False, sep='\t')


write_patient_dict_to_file(output_file)