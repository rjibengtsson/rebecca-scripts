##################
# Import modules #
##################

import sys, os
import pandas as pd
from pathlib import Path
import csv
import argparse

####################
# Take input files #
####################

parser = argparse.ArgumentParser(description='Format fragment length file to compare step2.sh output')
parser.add_argument('acc', help='Accession')
args = parser.parse_args()

sample_acc = args.acc


#############
# Functions #
#############

def format_file(sample_acc):
    # Store fragment length in a dictionary
    len_dict = {}
    acc = sample_acc
    pathlist = Path("./").glob('*.strain.length')
    for path in pathlist:
        path_in_str = str(path)
        df = pd.read_csv(path_in_str, names=['sample', 'length'], sep="\t")
        len = df.loc[df['sample'] == acc, 'length'].to_string(index=False)
        name = os.path.basename(path_in_str).split(".")[4:6]
        len_dict['_'.join(name)] = len

    # Create empty table
    table = pd.DataFrame(columns=['length', 'gte', 'pid'])

    # Format dictionary into table
    for key, value in len_dict.items():
        len = value
        gte = key.split("_")[0]
        pid = key.split("_")[1]
        table = table.append({'length': len, 'gte': gte[3:], 'pid': pid[3:]}, ignore_index=True)

    # Save table and execute R script
    table.to_csv(acc + '_fragLen.txt', sep='\t', index=False)

    # Generat R plot
    R_cmd = "Rscript frag_len_curve.R " + acc + '_fragLen.txt ' + acc + '_fragLen.png'
    os.system(R_cmd)


########
# Main #
########

format_file(sample_acc)