'''
This script will take a distance matrix and convert it into a pairwise distance list.
'''

import sys
import numpy as np
import pandas as pd
import io

in_handle = open(sys.argv[1]).read().splitlines()

array_lst = []
#
# for line in in_handle:
#     new_line = line.replace('\t', ' ')
#     array_lst.append(new_line)

txt_dist_mat = "\n".join(array_lst)
df = pd.read_fwf(io.StringIO(txt_dist_mat), index_col=0)
print(in_handle)

# array_lst=[]
#
# for line in in_handle:
#     new_line = str(line.replace('\t', ' '))
#     array_lst.append(new_line)
#
# txt_dist_mat = "\n".join(array_lst)
# print(array_lst)

#
# print(df)
#print(matrix)


#euc_mat = df.stack().reset_index().rename(columns={'level_0': 'Obj1', 'level_1': 'Obj2', 0: s