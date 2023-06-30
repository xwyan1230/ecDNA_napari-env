import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import numpy as np
import matplotlib as mpl
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 4

# samples
sample = 'DM_3_49pos'
figure_name = 'DM'
df = pd.read_csv('%s%s_summary.txt' % (data_dir, sample), na_values=['.'], sep='\t')
df['total_int_MYC'] = df['area_nuclear_IF'] * df['mean_int_MYC']
df['total_int_ecDNA'] = df['mean_int_ecDNA'] * df['total_area_ecDNA']
df['total_int_DNAFISH'] = df['mean_int_DNAFISH'] * df['area_nuclear']
df['total_int_hoechst'] = df['area_nuclear_IF'] * df['mean_int_hoechst']
df['ln_total_area_ecDNA'] = np.log(df['total_area_ecDNA'])
df['ln_total_int_MYC'] = np.log(df['total_int_MYC'])
df['ln_total_int_ecDNA'] = np.log(df['total_int_ecDNA'])
df['ln_total_int_DNAFISH'] = np.log(df['total_int_DNAFISH'])
df['ln_total_int_hoechst'] = np.log(df['total_int_hoechst'])

hue_order = ['unsync', 'G1', 'S', 'G2']

df_sample = df[df['cellcycle'].isin(hue_order)].copy().reset_index(drop=True)
df = df_sample
print(len(df))

plt.subplots(figsize=(12, 9))
x = 'ln_total_int_DNAFISH'
y = 'ln_total_int_hoechst'
sns.scatterplot(data=df, x=x, y=y, s=5)

"""hue_order1 = ['17.375', '17.625', '17.875', '18.125', '18.375', '18.625', '18.875']
cmap = mpl.cm.get_cmap('Spectral')
l = np.arange(0, 1, 1/len(hue_order1))
line_colors = [cmap(i) for i in l]
for k in range(len(hue_order1)):
    cutoff = float(hue_order1[k])
    data = df[(df['ln_total_int_MYC'] >= cutoff - 0.125) & (df['ln_total_int_MYC'] < cutoff + 0.125)].copy().reset_index(drop=True)
    sns.scatterplot(data=data, x=x, y=y, s=10, color=line_colors[k])
plt.savefig('%s/%s_%s_vs_%s_color-MYC.pdf' % (output_dir, figure_name, x, y))
plt.show()"""

hue_order1 = ['G1', 'S', 'G2']
line_colors = [(220/255, 20/255, 60/255), (154/255, 205/255, 50/255), (255/255, 127/255, 80/255)]
for k in range(len(hue_order1)):
    data = df[df['cellcycle'] == hue_order1[k]].copy().reset_index(drop=True)
    sns.scatterplot(data=data, x=x, y=y, s=10, color=line_colors[k])
plt.savefig('%s/%s_%s_vs_%s_color-cellcycle.pdf' % (output_dir, figure_name, x, y))
plt.show()

"""hue_order1 = ['G1', 'S', 'G2']
line_colors = [(220/255, 20/255, 60/255), (154/255, 205/255, 50/255), (255/255, 127/255, 80/255)]
for k in range(len(hue_order1)):
    data = df[df['cellcycle'] == hue_order1[k]].copy().reset_index(drop=True)
    sns.scatterplot(data=df, x=x, y=y, s=10, color='black', alpha=0.5)
    sns.scatterplot(data=data, x=x, y=y, s=10, color=line_colors[k])
    plt.savefig('%s/%s_%s_vs_%s_color-cellcycle_%s.pdf' % (output_dir, figure_name, x, y, hue_order1[k]))
    plt.show()"""
