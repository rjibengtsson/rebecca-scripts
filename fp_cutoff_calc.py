import sys
import numpy as np
import pandas as pd

# Read FP stats table into dataframe
df = pd.read_csv(sys.argv[1], header=0)
# print(df.shape)

filt_dict = {}
filt_dict['FPves / Kbp'] = []
filt_dict['Total FPves rmv'] = []
filt_dict['% FPves rmv'] = []
filt_dict['No. Codex blocks rmv'] = []
filt_dict['size of codex rmv (bp)'] = []
filt_dict['% of Codex rmv'] = []

# Define threshold range
threshold_range = np.round(np.arange(0.5,2.1,0.1), 2).tolist()

# Iterate through range
for i in threshold_range:

    # Filter column by FPves per Kb
    freq = df[df['FPve per Kb']>=i]

    # Append info to dictionary
    filt_dict['FPves / Kbp'].append(i)
    filt_dict['Total FPves rmv'].append(freq['N of Fpves'].sum())
    filt_dict['% FPves rmv'].append((freq['N of Fpves'].sum()/df['N of Fpves'].sum())*100)
    filt_dict['No. Codex blocks rmv'].append(freq.shape[0])
    filt_dict['size of codex rmv (bp)'].append(freq['Length Block'].sum())
    filt_dict['% of Codex rmv'].append((freq['Length Block'].sum())/6385295*100)

df2 = pd.DataFrame.from_dict(filt_dict)
df2.to_csv(sys.argv[2], index=False)







