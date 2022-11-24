import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

out_handle = open(sys.argv[2], 'w')

for record in SeqIO.parse(sys.argv[1], 'fasta'):
    record_header = record.id
    if record_header.startswith('03123017924237'):
        header_split = record_header.split('-')
        new_header = str(header_split[0] + "-" + header_split[1] + "-" + header_split[2])
        record.id = new_header
        record = SeqRecord(Seq(record.seq), id=record.id, name=record.id, description='', )
        SeqIO.write(record, out_handle, "fasta")
    else:
        record = SeqRecord(Seq(record.seq), id=record.id, name=record.id, description='', )
        SeqIO.write(record, out_handle, "fasta")
