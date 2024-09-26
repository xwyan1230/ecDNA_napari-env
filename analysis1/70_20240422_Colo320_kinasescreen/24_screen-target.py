import numpy as np
import shared.dataframe as dat
import pandas as pd
import collections
import os

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"

target_df = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % master_folder, na_values=['.'])
targets = [str(target_df['TargetGene'][i]).split(', ') for i in range(len(target_df))]
targets_temp = []
for i in range(len(targets)):
    for j in range(len(targets[i])):
        targets_temp.append(targets[i][j])
targets_set = set(targets_temp)
targets_duplicated = [x for x, count in collections.Counter(targets_temp).items() if count > 1]
print(len(targets_set))
print(targets_set)
print(len(targets_duplicated))
print(targets_duplicated)