'''
This script will split multi-fasta file into individual fasta files
'''

import sys,os
from Bio import SeqIO


counter = 0

fasta_name = sys.argv[1]

for record in SeqIO.parse(fasta_name, "fasta"):

    # set counter
    counter+=1

    # retrieve file name
    file_name = (os.path.splitext(fasta_name)[0] + '_' + str(counter) + '.fna')

    # write output
    SeqIO.write(record, file_name, "fasta")

