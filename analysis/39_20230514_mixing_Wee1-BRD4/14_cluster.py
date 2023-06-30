import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
from shared.sinaplot import sinaplot

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230514_analysis_mixing_Wee1-BRD4/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

exp = '20230512_mixing_Wee1-BRD4_6hr'
sample = '10_BRD4-GFP-6hr_Ctrl-mCh-6hr'
GFP_sample = 'BRD4_6hr'
mCherry_sample = 'Ctrl'
hue_order = [mCherry_sample, GFP_sample]

df = pd.read_csv('%s%s/%s_summary.txt' % (data_dir1, sample, sample), na_values=['.'], sep='\t')
sample_lst = []
for i in range(len(df)):
    if df['group'][i] == 'GFP':
        sample_lst.append(GFP_sample)
    elif df['group'][i] == 'mCherry':
        sample_lst.append(mCherry_sample)
    else:
        sample_lst.append('NA')
df['sample'] = sample_lst
df['r'] = np.sqrt(df['area_nuclear']/math.pi)
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA']/math.pi)
df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt']/df['r']
df['dis_to_hub_area_normalized'] = df['dis_to_hub_area']/df['r']
df_sort = df[df['total_area_ecDNA_sqrt_normalized'] > 0.15].copy().reset_index(drop=True)
df = df_sort

df_GFP = df[df['group'] == 'GFP'].copy().reset_index(drop=True)
df_mCherry = df[df['group'] == 'mCherry'].copy().reset_index(drop=True)

print(len(df))
print(len(df_GFP))
print(len(df_mCherry))

fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
# sns.violinplot(data=df_sort, x='MYC_group1', y='dis_to_hub_area_normalized', order=hue_order2)
sinaplot(data=df, x='sample', y='dis_to_hub_area_normalized', order=hue_order, violin=False, scale='area')
plt.savefig('%s%s/dis_to_hub_area_normalized.pdf' % (output_dir, sample))
plt.show()

print(np.mean(df_GFP['dis_to_hub_area_normalized']))
print(np.mean(df_mCherry['dis_to_hub_area_normalized']))

print(np.std(df_GFP['dis_to_hub_area_normalized']))
print(np.std(df_mCherry['dis_to_hub_area_normalized']))

print(np.quantile(df_GFP['dis_to_hub_area_normalized'].tolist(), np.arange(0, 1.1, 0.25))[2])
print(np.quantile(df_mCherry['dis_to_hub_area_normalized'].tolist(), np.arange(0, 1.1, 0.25))[2])

print(np.mean(df_GFP['n_ecDNA']))
print(np.mean(df_mCherry['n_ecDNA']))
