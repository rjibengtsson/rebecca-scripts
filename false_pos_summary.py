'''
This script will take a vcf file containing false positives and produce a table containing various metrics
'''

##################
# Import modules #
##################

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='This script will produce various false positive metrics')
parser.add_argument('-vcf', '--vcf', help='VCF file containing false positives', required=True)
parser.add_argument('-outfile', '--outfile', help='Output csv file', required=True)
parser.add_argument('-strains', '--strains', help='Number of strains in codex', type=int, required=True)
parser.add_argument('-codexsize', '--codexsize', help='Codex size', type=int, required=True)

args = vars(parser.parse_args())

VCFFILE = args.get('vcf')
OUTFILE = args.get('outfile')
STRAINnum = args.get('strains')
CODEXsize = args.get('codexsize')

#############
# Functions #
#############

def get_metrics(VCFFILE, OUTFILE, STRAINnum, CODEXsize):

    # Grab data from VCF file
    df = pd.read_csv(VCFFILE, header=None, sep="\t")

    final_dict = {}

    s = df.groupby([0]).size().to_frame()

    for index, row in s.iterrows():
        contig_header = index
        final_dict[contig_header] = []
        fp_num = row[0]
        num_strains = str(contig_header).split("-")[-2].replace('S', '')
        block_len = str(contig_header).split("-")[-1].replace('L', '')
        final_dict[contig_header].extend([fp_num, num_strains, block_len])

    new_dict = {}
    new_dict['Codex block name'] = []
    new_dict['Number of FPves'] = []
    new_dict['No. Strains in block'] = []
    new_dict['Block len'] = []
    new_dict['FPves per Kb'] = []
    new_dict['Impact'] = []

    for k, v in final_dict.items():
        fp_num, num_strains, block_len = v
        fp_per_kb = int(fp_num) / int(block_len) * 1000
        impact = (int(block_len) * int(num_strains)) / (STRAINnum * CODEXsize)
        new_dict['Codex block name'].append(k)
        new_dict['Number of FPves'].append(fp_num)
        new_dict['No. Strains in block'].append(num_strains)
        new_dict['Block len'].append(block_len)
        new_dict['FPves per Kb'].append(fp_per_kb)
        new_dict['Impact'].append(impact)
        # data = [k, fp_num, num_strains, block_len, fp_per_kb, impact]

    df2 = pd.DataFrame.from_dict(new_dict)

    # write dataframe to csv
    df2.to_csv(OUTFILE, index=False)

########
# Main #
########

get_metrics(VCFFILE, OUTFILE, STRAINnum, CODEXsize)