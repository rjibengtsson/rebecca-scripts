'''
This script grab kraken results from the customer pipeline QC.log file

Usage: python3 grab_kraken_results_from_QClog.py -list samples.txt -out kraken_QC.txt
'''

##################
# Import modules #
##################

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Grab kraken results from QC.log file')
parser.add_argument('-list','--list', help='List containing sample names', required=True)
parser.add_argument('-out','--out', help='Out file containing concatenated kraken results', required=True)
args = parser.parse_args()

sample_lst = args.list
out_file = args.out

def format_kraken_dict(sample_lst):

    samples = open(sample_lst, 'r').read().splitlines()

    kraken_dict = {}

    for sample in samples:
        file = f"./{sample}/QC.log"
        df = pd.read_csv(file, header=None, sep="\t")
        Kraken_Percent_Species = df.loc[1,1]
        Kraken_Percent_Genus = df.loc[2,1]
        kraken_dict[sample] = [Kraken_Percent_Species, Kraken_Percent_Genus]

    return kraken_dict

def write_dict_to_table(sample_lst, out_file):

    kraken_dict = format_kraken_dict(sample_lst)

    with open(out_file, 'w') as out:
        title = "Sample" + "\t" + "Kraken_Percent_Species" + "\t" + "Kraken_Percent_Genus" + "\n"
        out.write(title)
        for key, value in kraken_dict.items():
            line = key + "\t" + value[0] + "\t" + value[1] + "\n"
            out.write(line)



write_dict_to_table(sample_lst,out_file)