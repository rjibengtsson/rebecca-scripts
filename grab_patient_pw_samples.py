'''
Usage: python3 ~/scripts/python/grab_patient_pw_samples.py patient_pw_samples.txt patient_pw_samples.csv
'''

import sys
import pandas as pd

patient_id_lst = list(open(sys.argv[1], "r").read().splitlines()) # read patient id as list

df = pd.read_csv(sys.argv[2], header=0)

patient_dict = {}

# create dictionary
for patient in patient_id_lst:
    for i, r in df.iterrows():
        patient_id = r['Patient ID']
        acc_1 = r['AccessionID1']
        acc_2 = r['AccessionID2']
        if patient == patient_id:
            list_of_acc = patient_dict.get(patient, [])
            list_of_acc.append(acc_1)
            list_of_acc.append(acc_2)
            patient_dict[patient] = list_of_acc

# write into separate lists
for k, v in patient_dict.items():
    v = set(v)
    with open(k + '_vcf_list.txt', 'w') as file:
        file.write("\n".join(str('/scratch/rjibengtsson/ludden21/vcf_lowerMQ/' + acc + '.vcf.gz') for acc in v))