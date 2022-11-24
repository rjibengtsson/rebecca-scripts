'''
Author: Rebecca

Description: This script will take a blacklist bed file and a merged vcf file, remove variants from the vcf file in the
blacklist region and generate a new filtered vcf file, create a dendrogram based on the filtered vcf and generate stats on
how many sites are removed.

Usage: python3 blacklist_filtering.py -vcf ST7095_cluster1_merged_onlyP1.vcf.gz -bed 2_ec.bed -prefix doubleMapping
'''


##################
# Import modules #
##################

import textwrap
import argparse
import subprocess
import pathlib
import pandas as pd
import os, sys


####################
# Take input files #
####################

def parse_args():
    # set up command line arguments
    parser = argparse.ArgumentParser(
        prog = 'blacklist_filtering.py',
        usage='%(prog)s -vcf VCF -bed BED -prefix PREFIX',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        Perform blacklist filtering of vcf file and generate a dendrogram using RapidNJ
        '''))

    parser.add_argument('-vcf','--vcf', help='VCF file to be filtered', required=True)
    parser.add_argument('-bed','--bed', help='Black list bed file', required=True)
    parser.add_argument('-prefix','--prefix', help='Prefix of your output files', required=True)

    try:
        args = parser.parse_args()

        if not os.path.exists(args.vcf):
            sys.stderr.write(f"ERROR:Please check input file path for {args.vcf}\n")
            sys.exit(1)

        if not os.path.exists(args.bed):
            sys.stderr.write(f"ERROR:Please check input file path for {args.bed}\n")
            sys.exit(1)


    except Exception as e:
        sys.stderr.write(str(e))
        sys.exit(1)

    return args


#############
# Functions #
#############


def mask_vcf_file(args):

    cmd1 = f"mkdir {path}"
    subprocess.run([cmd1], shell=True)

    cmd2 = f"bedtools intersect -header -v -a {args.vcf} -b {args.bed} > {path}/{vcf_name}_{args.prefix}.vcf"
    subprocess.run([cmd2], shell=True)

    cmd3 = f"bgzip {path}/{vcf_name}_{args.prefix}.vcf"
    subprocess.run([cmd3], shell=True)


def run_plink2(args):

    # generate pgen
    cmd1 = f"plink2 --make-pgen --out {path}/plink2_output --vcf {path}/{vcf_name}_{args.prefix}.vcf.gz " \
           f"--snps-only just-acgt --set-missing-var-ids @:#\$r,\$a --allow-extra-chr"
    subprocess.run([cmd1], shell=True)

    # generate the kinship table counts
    cmd2 = f"plink2 --pfile {path}/plink2_output --snps-only just-acgt " \
           f"--set-missing-var-ids @:#\$r,\$a --allow-extra-chr " \
           f"--make-king-table counts cols=id,nsnp,ibs0 --threads 4"
    subprocess.run([cmd2], shell=True)

    # move kinship table into output directory
    cmd3 = f"mv plink2.kin0 plink2.log {path}/"
    subprocess.run([cmd3], shell=True)


def convert_plink2_distance_to_matrix():

    # convert plink2 distance to matrix
    cmd1 = f"/home/lmontemayor/Scripts/plink2matrix {path}/plink2.kin0"
    subprocess.run([cmd1], shell=True)

    # # move tsv files
    # cmd2 = f"mv final_normalised_matrix_RapidNJ.tsv kin0_matrix.tsv sdiff_matrix.tsv {path}/"
    # subprocess.run([cmd2], shell=True)

    return f'{path}/plink2.kin0_final_normalized_matrix.phylip'


def round_to_integers():

    rapidNJ_file = convert_plink2_distance_to_matrix()
    input_file = open(rapidNJ_file, 'r').readlines()
    out_file = f"{path}/final_normalised_matrix_int_RapidNJ.tsv"

    with open(f"{path}/tmp.tsv", 'w') as tmp_file:
        for line in input_file[1:]:
            tmp_file.write(line)

    df = pd.read_csv(f"{path}/tmp.tsv", sep="\t", header=None)

    column_index_len = len(df[0])
    df2 = df.set_index(df.columns[0])
    df3 = df2.round(0)

    df3.to_csv(f"{path}/tmp.tsv", sep='\t')

    with open(f"{path}/tmp.tsv", 'r') as data_file:
        lines = data_file.readlines()
        lines[0] = str(column_index_len) + "\n"
        with open(out_file, 'w') as out_data:
            for line in lines:  # write updated lines
                out_data.write(line)

    # remove tmp file
    cmd = f'rm {path}/tmp.tsv'
    subprocess.run([cmd], shell=True)

    return out_file


def generate_dendrogram(args):

    # input_file = round_to_integers()
    input_file = convert_plink2_distance_to_matrix()

    cmd = f"/home/lmontemayor/bin/rapidNJ/bin/rapidnj {input_file} -i pd -n -x {path}/{args.prefix}_filtered.nwk"
    subprocess.run([cmd], shell=True)

    return f"{path}/{args.prefix}_filtered.nwk"


def generate_stats(args):

    cmd = f"bash ~/software/scripts/count_black_listing_sites.sh -b {args.bed} " \
          f"-v {args.vcf} -f {path}/{vcf_name}_{args.prefix}.vcf.gz " \
          f"-o {path}/number_of_masked_sites.txt"

    subprocess.run([cmd], shell=True)

    return f"{path}/number_of_masked_sites.txt"

########
# Main #
########

if __name__ == "__main__":
    args = parse_args()

    # Define vcf name and path
    vcf_name = str(args.vcf).split(".")[0]
    path = pathlib.Path.cwd().joinpath(f'{args.prefix}')

    # Call functions
    mask_vcf_file(args)
    run_plink2(args)
    convert_plink2_distance_to_matrix()
    round_to_integers()
    generate_dendrogram(args)
    generate_stats(args)