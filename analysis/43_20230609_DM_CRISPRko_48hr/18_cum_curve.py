import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
from shared.sinaplot import sinaplot
import shared.dataframe as dat
import seaborn as sns
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 4

# samples
sample = 'C2'
hue_order = ['GFP', 'mCherry']

df = pd.read_csv('%s/%s_n4_simplified.txt' % (data_dir, sample), na_values=['.'], sep='\t')
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
feature = ['percentage_area_curve_ecDNA', 'cum_area_ind_ecDNA_filled']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]
df_GFP = df[df['group'] == 'GFP'].copy().reset_index(drop=True)
df_mCherry = df[df['group'] == 'mCherry'].copy().reset_index(drop=True)

# cumulative curve
print("Plotting cumulative curve...")

for f in feature:
    x_label = 'number of ecDNA hub'
    df1 = df_GFP
    df2 = df_mCherry
    number_nuclear1 = len(df1)
    number_nuclear2 = len(df2)

    mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(dat.list_fill_with_last_num(df1[f].tolist()))
    mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(dat.list_fill_with_last_num(df2[f].tolist()))

    plt.subplots(figsize=(6, 4))
    for i in range(len(df1)):
        plt.plot(df1[f][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
    for i in range(len(df2)):
        plt.plot(df2[f][i], alpha=0.05, color=[line_colors[1][j]+0.05 for j in range(len(line_colors[1]))])
    plt.plot(mean_curve1, color=line_colors[0], label='%s, n=%s' % (hue_order[0], number_nuclear1))
    plt.plot(mean_curve2, color=line_colors[1], label='%s, n=%s' % (hue_order[1], number_nuclear2))
    plt.plot(ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.xlabel(x_label)
    plt.ylabel(f)
    plt.legend()
    plt.savefig('%s/%s_%s.pdf' % (output_dir, sample, f))
    plt.show()