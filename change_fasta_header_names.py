import sys
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# read header files
df = pd.read_csv(sys.argv[1], header=None, sep="\t")

out_handle = open(sys.argv[3], 'w')

for index, row in df.iterrows():
    for record in SeqIO.parse(sys.argv[2], 'fasta'):
        if record.id == row[0]:
            record.id = row[3]
            record = SeqRecord(Seq(record.seq), id=record.id, name=record.id, description='', )
            SeqIO.write(record, out_handle, "fasta")