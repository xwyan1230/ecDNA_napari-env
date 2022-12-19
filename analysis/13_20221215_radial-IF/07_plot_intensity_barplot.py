import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import matplotlib as mpl
import numpy as np
from matplotlib import cm
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221215_analysis_radial-IF/20221117_immunoFISH_acid/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'H3K4me3'
n = 8
df = pd.read_csv(("%s%s_radial_n%s.txt" % (data_dir, sample, n)), na_values=['.'], sep='\t')

line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25),  (0.30, 0.70, 0.70)]

feature = ['radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_IF',
           'radial_curve_DNAFISH_normalized', 'radial_curve_IF_normalized',
           'nuclear_int', 'DNAFISH_int', 'IF_int', 'ecDNA_label',
           'int_r_to_edge', 'int_r_to_center', 'int_relative_r']

for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

x = np.arange(0.025, 1, 0.05)
x_str = ['%.3f' % elem for elem in x]
df_bar = pd.DataFrame()

for k in x:
    IF_ecDNA_neg_lst = []
    IF_ecDNA_pos_lst = []
    for i in range(len(df)):
        temp = [df['IF_int'][i][j] for j, e in enumerate(df['int_relative_r'][i]) if (e > (k-0.025)) & (e <= (k+0.025))]
        ecDNA = [df['ecDNA_label'][i][j] for j, e in enumerate(df['int_relative_r'][i]) if (e > (k-0.025)) & (e <= (k+0.025))]
        IF_ecDNA_neg_lst = IF_ecDNA_neg_lst + [temp[j] for j, e in enumerate(ecDNA) if e == 0]
        IF_ecDNA_pos_lst = IF_ecDNA_pos_lst + [temp[j] for j, e in enumerate(ecDNA) if e > 0]
    print(len(IF_ecDNA_neg_lst))
    print(len(IF_ecDNA_pos_lst))
    df_bar_temp = pd.DataFrame({'sample': ['ecDNA_pos'] * len(IF_ecDNA_pos_lst) + ['ecDNA_neg'] * len(IF_ecDNA_neg_lst),
                                'IF_int': IF_ecDNA_pos_lst + IF_ecDNA_neg_lst,
                                'relative_r': [k] * (len(IF_ecDNA_pos_lst)+len(IF_ecDNA_neg_lst))})
    df_bar = pd.concat([df_bar, df_bar_temp])
df_bar.to_csv('%s%s_int_n%s.txt' % (output_dir, sample, n), index=False, sep='\t')

# plot bar plot
sns.set_palette(sns.color_palette(line_colors))
hue_order = ['ecDNA_neg', 'ecDNA_pos']
fig, ax = plt.subplots(figsize=(12, 9))
fig.subplots_adjust(left=0.2)
ax = sns.barplot(data=df_bar, x='relative_r', y='IF_int', estimator=np.mean, hue='sample', hue_order=hue_order)
plt.savefig('%s/barplot_relative_r_%s.pdf' % (output_dir, sample))
plt.close()