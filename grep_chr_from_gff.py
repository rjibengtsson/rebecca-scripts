import os, sys

in_handle = open(sys.argv[1]).read().splitlines()

for line in in_handle:
    refseq = line.split("\t")[0]
    acc = line.split("\t")[1]
    cmd = "grep " + '"' + acc + '" ' + refseq + ".gff > " + refseq + ".chr.gff "
    os.system(cmd)