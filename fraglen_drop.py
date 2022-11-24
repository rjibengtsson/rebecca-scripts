'''
This script takes the outputs of step2.sh and visualise % fragment length drop off at different pid threshold
'''

##################
# Import modules #
##################

import sys, os
import pandas as pd
from glob import glob
import argparse

# sort files by pid values
file_lst = []
gte = 21
path_to_search = f"/home/rjibengtsson/codex_ecoli/phase2/ecor/phase2/all.exclISE.missing.fasta.gte{gte}.pid*.strain.length"
target_files = glob(path_to_search)

for path in target_files:
    path_in_str = str(path)
    name = os.path.basename(path_in_str).split(".")[5].replace('pid', '')
    file_lst.append(name)

file_lst_srt = sorted(file_lst)

# iterate through files and store length as dictionary
length_dic = {}

for file in file_lst_srt:
    path = f"/home/rjibengtsson/codex_ecoli/phase2/ecor/phase2/all.exclISE.missing.fasta.gte{gte}.pid{file}.strain.length"
    df = pd.read_csv(path, names=['sample', 'length'], sep="\t")
    for index, row in df.iterrows():
        if row['sample'] in length_dic.keys():
            length_dic[row['sample']].append(row['length'])
        else:
            length_dic[row['sample']] = [row['length']]

# write new dataframe
new_df = pd.DataFrame.from_dict(length_dic, orient='index')
new_df = new_df.assign(diff_1 = new_df[0] - new_df[1])
new_df = new_df.assign(diff_2 = new_df[1] - new_df[2])

new_df.to_csv(sys.argv[1])