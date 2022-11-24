import sys
from Bio import SeqIO

fasta_file = sys.argv[1]
fasta_file_name =  str(fasta_file).split('.')[0]

for record in SeqIO.parse(fasta_file, 'fasta'):
    SeqIO.write(record, fasta_file_name + '.fasta', "fasta")