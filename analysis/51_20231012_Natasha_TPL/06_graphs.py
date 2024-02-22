import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from shared.sinaplot import sinaplot

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230724_analysis_Natasha_colcemid/TPL/mCh-DMSO_GFP-TPL/"
output_dir = "%s06_analysis/" % master_folder

sample = 'mCh-DMSO_GFP-TPL'

df = pd.read_csv('%s/06_analysis/%s_1_n2.txt' % (master_folder, sample), na_values=['.'], sep='\t')
df1 = pd.read_csv('%s/06_analysis/%s_3_n2.txt' % (master_folder, sample), na_values=['.'], sep='\t')
df2 = pd.read_csv('%s/06_analysis/%s_4_n2.txt' % (master_folder, sample), na_values=['.'], sep='\t')
df = pd.concat([df, df1, df2], axis=0).reset_index(drop=True)

hue_order = ['mCherry', 'GFP']
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
sns.set_palette(sns.color_palette(line_colors))

df['ln_mean_int_green'] = np.log(df['mean_int_green'])
df['ln_mean_int_red'] = np.log(df['mean_int_red'])
df['r'] = np.sqrt(df['area_nuclear']/math.pi)
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA']/math.pi)
df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt']/df['r']
df['dis_to_hub_area_normalized'] = df['dis_to_hub_area']/df['r']
df_sort = df[df['total_area_ecDNA_sqrt_normalized'] > 0.1].copy().reset_index(drop=True)
df = df_sort

sample_lst = []
for i in range(len(df)):
    if (df['ln_mean_int_red'][i] < 7) & (df['ln_mean_int_green'][i] > 8.75):
        sample_lst.append('GFP')
    elif (df['ln_mean_int_red'][i] > 8.75) & (df['ln_mean_int_green'][i] < 7):
        sample_lst.append('mCherry')
    else:
        sample_lst.append('NA')
df['group'] = sample_lst

df_GFP = df[(df['ln_mean_int_green'] > 8.75) & (df['ln_mean_int_red'] < 7)].copy().reset_index(drop=True)
df_mCherry = df[(df['ln_mean_int_green'] < 8.75) & (df['ln_mean_int_red'] > 7)].copy().reset_index(drop=True)


print(len(df))
print(len(df_GFP))
print(len(df_mCherry))
print(np.mean(df_GFP['dis_to_hub_area_normalized']))
print(np.mean(df_mCherry['dis_to_hub_area_normalized']))

fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
# sns.violinplot(data=df_sort, x='MYC_group1', y='dis_to_hub_area_normalized', order=hue_order2)
sinaplot(data=df_sort, x='group', y='dis_to_hub_area_normalized', order=hue_order, violin=False, scale='area')
plt.savefig('%s/%s_dis_to_hub_area_normalized.pdf' % (output_dir, sample))
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_GFP, y='dis_to_hub_area_normalized', x='total_area_ecDNA_sqrt_normalized', color=line_colors[0], s=10, alpha=1)
sns.scatterplot(data=df_mCherry, y='dis_to_hub_area_normalized', x='total_area_ecDNA_sqrt_normalized', color=line_colors[1], s=10, alpha=1)
plt.legend()
plt.savefig('%s/%s_dis_to_hub_area_normalized_scatter.pdf' % (output_dir, sample))
plt.show()
