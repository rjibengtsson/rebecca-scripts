import sys
import pandas as pd

# read file as pandas dataframe
df = pd.read_csv(sys.argv[1], header=None, sep="\t")

# iterate through dataframe
for index, row in df.iterrows():
    if str(row["#CHROM"]).startswith('0314'):
        name_split = str(row["#CHROM"]).split('-')
        new_name = str(name_split[0] + "-" + name_split[1] + "-" + name_split[2] + "-" + name_split[3])
        df.at[index, "#CHROM"] = new_name
    else:
        continue

df.to_csv(sys.argv[2], sep='\t',index=False)
