##################
# Import modules #
##################

import sys
from Bio import SeqIO
import re

locus_lst = list(open(sys.argv[1], 'r').read().splitlines())
out_handle = open(sys.argv[3], 'w')

for locus_id in locus_lst:
    locus_tag = locus_id
    pattern = re.compile('.*({}).*'.format(locus_tag))
    for record in SeqIO.parse(sys.argv[2], 'fasta'):
        record_name = str(record.description)
        if pattern.match(record_name):
            # record.id = record_name
            # record.description = record_name
            SeqIO.write(record, out_handle, 'fasta')
        else:
            pass