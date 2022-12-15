import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221215_analysis_MYC-IF/20221117_immunoFISH_acid/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'MYC-old'
df = pd.read_csv(("%s%s_ecDNA.txt" % (data_dir, sample)), na_values=['.'], sep='\t')

hist_colors = [(0.90, 0.90, 0.90), (0.95, 0.50, 0.50),  (0.50, 0.90, 0.90)]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25),  (0.30, 0.70, 0.70)]

df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA'])
df['dis_to_hub_area_normalized'] = df['dis_to_hub_area']/np.sqrt(df['area_nuclear'])

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, x='total_area_ecDNA', y='dis_to_hub_area', c=df['MYC_mean'], vmax=40000)
plt.savefig('%s/dis_to_hub_area_vs_total_area_%s.pdf' % (output_dir, sample))
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, x='total_area_ecDNA_sqrt', y='dis_to_hub_area', c=df['MYC_mean'], vmax=40000)
plt.savefig('%s/dis_to_hub_area_vs_total_area_sqrt_%s.pdf' % (output_dir, sample))
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, x='total_area_ecDNA_sqrt', y='dis_to_hub_area_normalized', c=df['MYC_mean'], vmax=40000)
plt.savefig('%s/dis_to_hub_area_normalized_vs_total_area_sqrt_%s.pdf' % (output_dir, sample))
plt.show()