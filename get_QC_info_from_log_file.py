'''
Description: This script will grab various info from the QC.log customer pipeline output files, format it into a
tab-delimited file and draw frequency distribution plots based on various QC info

Usage: python get_QC_info_from_log_files.py -list samples.txt -out QC_info.txt
'''

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Grab info from QC.log file and draw frequency distribution plots')
parser.add_argument('-list','--list', help='List containing sample names', required=True)
parser.add_argument('-out','--out', help='Tab-delimited out file containing concatenated QC results', required=True)
args = parser.parse_args()

sample_lst = args.list
out_file = args.out


def format_QC_dict(sample_lst):

    samples = open(sample_lst, 'r').read().splitlines()

    QC_dict = {}

    for sample in samples:
        file = f"./{sample}.QC.log"
        df = pd.read_csv(file, header=None, sep="\t")
        Kraken_Percent_Species = df.loc[1,1]
        Kraken_Percent_Genus = df.loc[2,1]
        Raw_Read_Count = df.loc[3,1]
        Raw_Base_Count = df.loc[4,1]
        Trimmed_Read_Count = df.loc[5,1]
        Trimmed_Base_Count = df.loc[6,1]
        QC_dict[sample] = [Kraken_Percent_Species, Kraken_Percent_Genus, Raw_Read_Count, Raw_Base_Count, Trimmed_Read_Count, Trimmed_Base_Count]

    return QC_dict

def write_dict_to_table(sample_lst, out_file):

    QC_dict = format_QC_dict(sample_lst)

    with open(out_file, 'w') as out:
        title = "Sample" + "\t" + "Kraken_Percent_Species" + "\t" + "Kraken_Percent_Genus" + "\t" + \
                "Raw_Read_Count" + "\t" + "Raw_Base_Count" + "\t" + "Trimmed_Read_Count" + "\t" + "Trimmed_Base_Count" + "\n"
        out.write(title)
        for key, value in QC_dict.items():
            line = key + "\t" + value[0] + "\t" + value[1] + "\t" + value[2] + "\t" + value[3] + "\t" + value[4] + "\t" + value[5] + "\n"
            out.write(line)

    return out_file

# def plot_distributions():
#
#     qc_out_file = write_dict_to_table(sample_lst, out_file)
#

write_dict_to_table(sample_lst, out_file)