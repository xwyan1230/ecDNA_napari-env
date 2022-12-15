import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221213_analysis_DNA-FISH_mixing_exp/20221121_H2B-series_POLR3D-BRD4-BRD1-DAPK2/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'H2B+POLR3D'
pos_threshold = 20000
neg_threshold = 12000
pos = 'ctrl'
neg = 'POLR3D KO'
hue_order = [pos, neg]

df = pd.read_csv(("%s%s_ecDNA.txt" % (data_dir, sample)), na_values=['.'], sep='\t')

hist_colors = [(0.90, 0.90, 0.90), (0.95, 0.50, 0.50),  (0.50, 0.90, 0.90)]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25),  (0.30, 0.70, 0.70)]

sample_lst = []
for i in range(len(df)):
    if df['mCherry_mean'].tolist()[i] < neg_threshold:
        sample_lst.append(neg)
    elif df['mCherry_mean'].tolist()[i] > pos_threshold:
        sample_lst.append(pos)
    else:
        sample_lst.append('NA')
df['sample'] = sample_lst

df_sort = df[df['sample'].isin([neg, pos])].copy().reset_index(drop=True)
df_pos = df[df['sample'].isin([pos])].copy().reset_index(drop=True)
df_neg = df[df['sample'].isin([neg])].copy().reset_index(drop=True)

df_sort['total_area_ecDNA_sqrt'] = np.sqrt(df_sort['total_area_ecDNA'])
df_sort['dis_to_hub_area_normalized'] = df_sort['dis_to_hub_area']/np.sqrt(df_sort['area_nuclear'])

"""sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_sort, x='total_area_ecDNA', y='dis_to_hub_area', hue='sample', hue_order=hue_order)
plt.savefig('%s/dis_to_hub_area_vs_total_area_%s.pdf' % (output_dir, sample))
plt.show()"""

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_sort, x='total_area_ecDNA_sqrt', y='dis_to_hub_area', hue='sample', hue_order=hue_order)
plt.savefig('%s/dis_to_hub_area_vs_total_area_sqrt_%s.pdf' % (output_dir, sample))
plt.show()
