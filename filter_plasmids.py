#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a script file to filter sequences of plasmids in downloaded genomes.
It looks for 'plasmid' in the sequence record and saves those sequences that do not have it.
"""

from Bio import SeqIO
import os
import argparse
import sys

fasta_file = sys.argv[1]
out_handle = sys.argv[2]
def remove_plasmids(fasta_file):
    for seq in SeqIO.parse(open(fasta_file), 'fasta'):
        if "plasmid" not in seq.description:
            out_fasta = open(out_handle, "w")
            SeqIO.write(seq, out_fasta, "fasta")

def main():
        filter_plasmids = remove_plasmids(fasta_file)


if __name__ == "__main__":
        main()