import sys
import pandas as pd

df = pd.read_csv(sys.argv[1], header=0, sep='\t')

PG_list = ['A', 'B1', 'B2', 'C', 'D', 'E', 'F', 'G', 'cladeI']

# create source type list
source_type_list = list(df['Source Type'].unique())

# create dictionary
final_dict ={}

# iterate through PGs
for PG in PG_list:
    PG_dict = {}
    df2 = df[df["Clermont Type (ClermonTyping)(Phylotypes)"] == PG]
    for source in source_type_list:
        source_type_count = (df2['Source Type'] == source).sum()
        PG_dict[source] = source_type_count
    final_dict[PG] = PG_dict

# print(final_dict['C'])

append_dict = {}

for key,value in final_dict.items():
    value_list = []
    for i in value:
        value_list.append(value[i])
    append_dict[key] = value_list

row_labels = source_type_list
df = pd.DataFrame(data=append_dict, index=row_labels)
# print(df)

df.to_csv(sys.argv[2], sep='\t')


