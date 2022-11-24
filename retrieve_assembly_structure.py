'''
Author: Rebecca

Description: This script will retrieve sequencing and assembly information from assembly structure report.

Usage: python3 retrieve_assembly_structure.py -dir ncbi-genomes-2022-11-22
'''


##################
# Import modules #
##################
import re
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
        prog = 'retrieve_assembly_structure.py',
        usage='%(prog)s -dir DIR',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        Retrieve sample sequencing and assembly metadata
        '''))

    parser.add_argument('-dir','--dir', help='Directory containing assembly structure reports', required=True)

    try:
        args = parser.parse_args()

        if not os.path.exists(args.dir):
            sys.stderr.write(f"ERROR:Please check path for {args.dir}\n")
            sys.exit(1)

    except Exception as e:
        sys.stderr.write(str(e))
        sys.exit(1)

    return args



column = ['Assembly name', 'Isolate', 'Infraspecific name', 'BioSample', 'BioProject', 'GenBank assembly accession',
          'RefSeq assembly accession', 'Date', 'Sequencing technology', 'Assembly method', 'Genome coverage']


def set_directory_path(args):

    directory_path = args.dir

    if directory_path.endswith('/'):
        new_directory_path = directory_path[:-1]
    else:
        new_directory_path = directory_path

    path = pathlib.Path.cwd().joinpath(f"{new_directory_path}/")

    return path


def report_to_list(in_file):

    list_of_lists = []

    with open(in_file, 'r') as file:
        for line in file:
            if line.startswith("# "):
                line = line.lstrip(line[0:2])
                line = line.strip()
                list_of_lists.append(line)

    return list_of_lists


def scan_list(in_list):

    search_lst = column

    assembly_dict ={}


    for key in search_lst:
        # print("key", key, assembly_dict)
        for item in in_list:
            # print("item", item)
            if key not in assembly_dict:
                if re.search(key, item):
                    assembly_dict[key] = [item.split(':')[1].strip()]
                    break
        else:
            assembly_dict[key] = ["NA"]

    return assembly_dict

def make_input_file_list():

    directory_path = set_directory_path(args)

    cmd = f"ls {directory_path}/*_assembly_report.txt > /{directory_path}/assembly_report_list.txt"

    subprocess.run([cmd], shell=True)

    return f"{directory_path}/assembly_report_list.txt"


def append_dict_to_df():

    file_list = open(make_input_file_list(), 'r').read().splitlines()

    df_assembly_dict = pd.DataFrame(columns=column)

    for file in file_list:
        list_of_lists = report_to_list(file)
        assembly_dict = scan_list(list_of_lists)
        df = pd.DataFrame(assembly_dict)
        df_assembly_dict = pd.concat([df_assembly_dict, df])

    directory_path = set_directory_path(args)

    df_assembly_dict.to_csv(f"{directory_path}/assembly_report.txt", sep="\t", index=False)



########
# Main #
########

if __name__ == "__main__":
    args = parse_args()

    # Call functions
    append_dict_to_df()