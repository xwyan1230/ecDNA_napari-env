import pandas as pd
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt
import seaborn as sns
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221017_periphery-localization_analysis/20220826_neuroblastoma/"
sample_lst = ['AC', 'AK', 'CB', 'EB', 'FR', 'GL', 'GO', 'IC', 'NE', 'V']
group = 'AC'
version = 1

df = pd.read_csv('%s%s_df_v%s.txt' % (master_folder, group, version), na_values=['.'], sep='\t')

plt.subplots(figsize=(12, 9))
plt.scatter(df[df['group'] == 'A']['radial_distance'], df[df['group'] == 'A']['k_cluster'], color='red')
plt.scatter(df[df['group'] == 'C']['radial_distance'], df[df['group'] == 'C']['k_cluster'], color='blue')
plt.xlabel('radial_distance')
plt.ylabel('k_cluster')
plt.legend()
# plt.savefig('%s/%s-sample-radial_distance-vs-k_cluster.pdf' % (save_folder, group))
plt.show()