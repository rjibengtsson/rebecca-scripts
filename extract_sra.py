'''
This script takes the File3_metadata_661K.txt metadata file from the 661k project paper
and selects a list of unique ST from each bio-project. The output of this script is a
tab-delimited file containing a list of BioSample id's with it's associated ena run accession numbers

Parameters:
    --
'''



import sys
import json
import re
import pandas as pd

df = pd.read_csv(sys.argv[1], header=0, delimiter='\t')

organism = "Escherichia coli"

new_df = df[df['scientific_name'].str.match(r'^' + organism + '.*') == True]

new_df.to_csv(sys.argv[2], sep='\t', index=False)


# with open(sys.argv[1], 'r') as json_file:
#     data = json.load(json_file)
#
# biosample_id_lst = list(open(sys.argv[2], 'r').read().splitlines())
# out_handle = open(sys.argv[3], 'w')
#
# for biosample_id in biosample_id_lst:
#     id = data[biosample_id]['ena_metadata'][0]
#     run_accession = id.get('run_accession')
#     line = biosample_id + "\t" + run_accession + "\n"
#     out_handle.write(line)
