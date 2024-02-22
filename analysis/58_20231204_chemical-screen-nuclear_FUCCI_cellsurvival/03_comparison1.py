import pandas as pd
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231204_analysis_FUCCI-screen/"
samples = ['DM_FUCCI_24hr', 'HSR_FUCCI_24hr']
ctrl = ['C3', 'C10', 'D6', 'F3', 'F10']

df_seq = pd.read_csv('%sseq.txt' % master_folder, na_values=['.'], sep='\t')

df = pd.read_csv('%s/%s/summary_survival.txt' % (master_folder, samples[0]), na_values=['.'], sep='\t')
df1 = pd.read_csv('%s/%s/summary_survival.txt' % (master_folder, samples[1]), na_values=['.'], sep='\t')

temp = []
for i in range(len(df_seq)):
    temp.append(df[df['sample'] == df_seq['well'][i]]['per_survival'].tolist()[0]/df1[df1['sample'] == df_seq['well'][i]]['per_survival'].tolist()[0])
df_feature = pd.DataFrame(temp)
df_feature.index = df_seq['compound'].tolist()
df_feature.columns = ['per_survival']

fig, ax = plt.subplots(figsize=(9, 14))
fig.subplots_adjust(left=0.3)
ax1 = sns.heatmap(df_feature, cbar=0, linewidths=2, vmax=2, vmin=0, square=True, cmap='coolwarm', annot=False, fmt='.2f') # annot=True
# plt.savefig('%s/24hr_survival_comparison.pdf' % master_folder)
plt.close()

df_sort = df_feature.sort_values(by='per_survival')
df_sort['chemical'] = df_sort.index
print(df_sort)
df_sort.to_csv('%s/24hr_per_survival_DM_vs_HSR.txt' % master_folder, index=False, sep='\t')
