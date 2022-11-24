import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

counter = 0
current_contig = None

output_handle = open(sys.argv[2], 'w')

with open(sys.argv[1]) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if current_contig != record.id:
            current_contig = record.id
            record.id = current_contig
            record = SeqRecord(Seq(record.seq),id=record.id, name=record.id, description='',)
            SeqIO.write(record, output_handle, "fasta")
            counter = 0
        else:
            record.id = current_contig + '_phage' + str(counter)
            record = SeqRecord(Seq(record.seq), id=record.id, name=record.id, description='', )
            SeqIO.write(record, output_handle, "fasta")
        counter += 1