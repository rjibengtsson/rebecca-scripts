import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

fasta_file = sys.argv[1]
fasta_file_name =  str(fasta_file).split('.')[0]

for record in SeqIO.parse(fasta_file, 'fasta'):
    seq = str(record.seq).replace('-', '')
    seq_A = seq.replace('a', 'A')
    seq_T = seq_A.replace('t', 'T')
    seq_G = seq_T.replace('g', 'G')
    new_seq = seq_G.replace('c', 'C')
    record = SeqRecord(Seq(new_seq), id=record.id, name=record.id, description='', )
    SeqIO.write(record, fasta_file_name + '.gap_stripped.fasta', "fasta")