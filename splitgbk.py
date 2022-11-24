from Bio import SeqIO
import sys

for rec in SeqIO.parse(sys.argv[1], "genbank"):
    print(rec)
   SeqIO.write([rec], open(rec.id + ".gbk", "w"), "genbank")
