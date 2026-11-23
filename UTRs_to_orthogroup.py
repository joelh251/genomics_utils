import pandas as pd
import gc

orthogroupsData = pd.read_csv("test_data/Orthogroups.tsv", sep="\t")

orthogroupDict = {}
for col in orthogroupsData.iloc[:, 1::]:
    orthogroupsData[col] = orthogroupsData[col].str.split(',')

    for i in range(len(orthogroupsData[col])):
        item = orthogroupsData[col].iloc[i] 
        if isinstance(item, list):
            for value in item:
                orthogroupDict[value.strip()] = orthogroupsData.iloc[i, 0]
        elif pd.notna(item):
            orthogroupDict[item.strip()] = orthogroupsData.iloc[i, 0]

del orthogroupsData
gc.collect()

UTRs = pd.read_csv("test_data/new_target_scan_UTRs.tsv", sep="\t", dtype={0: str}, keep_default_na=False, header=None)
UTRs = UTRs[UTRs.iloc[:, 0].str.strip() != 'NA']

#for i in range(len(UTRs.iloc[:, 0])):
#    item = UTRs.iloc[i, 0]
#    if pd.notna(item) and item in orthogroupDict:
#        UTRs.iloc[i, 0] = orthogroupDict[item]
UTRs.iloc[:, 0] = UTRs.iloc[:, 0].map(lambda x: orthogroupDict.get(x, x))
UTRs.to_csv("test_data/UTRs_test.tsv", sep="\t", index=False, header=None)