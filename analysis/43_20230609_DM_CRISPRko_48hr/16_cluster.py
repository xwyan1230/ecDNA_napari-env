import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from shared.sinaplot import sinaplot

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'C2'
hue_order = ['GFP', 'mCherry']
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
sns.set_palette(sns.color_palette(line_colors))

df = pd.read_csv('%s/%s/%s_n4_simplified.txt' % (data_dir1, sample, sample), na_values=['.'], sep='\t')
df['r'] = np.sqrt(df['area_nuclear']/math.pi)
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA']/math.pi)
df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt']/df['r']
df['dis_to_hub_area_normalized'] = df['dis_to_hub_area']/df['r']
df_sort = df[df['total_area_ecDNA_sqrt_normalized'] > 0.2].copy().reset_index(drop=True)
df = df_sort

df_GFP = df[df['group'] == 'GFP'].copy().reset_index(drop=True)
df_mCherry = df[df['group'] == 'mCherry'].copy().reset_index(drop=True)

print(len(df))
print(len(df_GFP))
print(len(df_mCherry))
print(np.mean(df_GFP['dis_to_hub_area_normalized']))
print(np.mean(df_mCherry['dis_to_hub_area_normalized']))

fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
# sns.violinplot(data=df_sort, x='MYC_group1', y='dis_to_hub_area_normalized', order=hue_order2)
sinaplot(data=df_sort, x='group', y='dis_to_hub_area_normalized', order=hue_order, violin=False, scale='area')
plt.savefig('%s/%s/%s_dis_to_hub_area_normalized.pdf' % (output_dir, sample, sample))
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_GFP, y='dis_to_hub_area_normalized', x='total_area_ecDNA_sqrt_normalized', color=line_colors[0], s=10, alpha=1)
sns.scatterplot(data=df_mCherry, y='dis_to_hub_area_normalized', x='total_area_ecDNA_sqrt_normalized', color=line_colors[1], s=10, alpha=1)
plt.legend()
plt.savefig('%s/%s/%s_dis_to_hub_area_normalized_scatter.pdf' % (output_dir, sample, sample))
plt.show()
