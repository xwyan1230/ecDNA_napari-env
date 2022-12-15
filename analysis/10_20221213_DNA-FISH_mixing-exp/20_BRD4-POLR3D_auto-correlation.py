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

hist_colors = [(0.90, 0.90, 0.90), (0.95, 0.50, 0.50),  (0.50, 0.90, 0.90)]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25),  (0.30, 0.70, 0.70)]

df = pd.read_csv(("%s%s_ecDNA.txt" % (data_dir, sample)), na_values=['.'], sep='\t')

feature = ['g']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA'])

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

# radial curve
print("Plotting g curve...")
x = np.arange(0, 101, 1)
x_label = 'r'

limit = 80

"""for s in hue_order:
    data = df[df['sample'] == s].copy().reset_index(drop=True)
    number_nuclear = len(data)

    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['g'].tolist())

    plt.subplots(figsize=(12, 9))
    for i in range(len(data)):
        plt.plot(x[:limit], data['g'][i][:limit], alpha=0.05, color=[line_colors[1][j] + 0.05 for j in range(len(line_colors[1]))])
    plt.plot(x[:limit], mean_curve3[:limit], color=line_colors[1], label='%s, n=%s' % ('ecDNA (normalized)', number_nuclear))
    plt.plot(x[:limit], ci_lower3[:limit], color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(x[:limit], ci_higher3[:limit], color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.axhline(y=1, color='black', linestyle='--')
    plt.xlabel(x_label)
    plt.ylim([0.5, 2])
    plt.ylabel('g')
    plt.legend()
    plt.savefig('%sg_%s_%s.pdf' % (output_dir, sample, s))
    plt.show()"""

plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['sample'] == hue_order[k]].copy().reset_index(drop=True)
    number_nuclear = len(data)

    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['g'].tolist())

    for i in range(len(data)):
        plt.plot(x[1:limit], data['g'][i][1:limit], alpha=0.05, color=[line_colors[k][j] + 0.05 for j in range(len(line_colors[k]))])
    plt.plot(x[1:limit], mean_curve3[1:limit], color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x[1:limit], ci_lower3[1:limit], color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x[1:limit], ci_higher3[1:limit], color=line_colors[k], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylabel('g')
plt.ylim([0.7, 1.6])
plt.legend()
plt.savefig('%sg_%s.pdf' % (output_dir, sample))
plt.show()

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_sort, x='total_area_ecDNA_sqrt', y='g_value', hue='sample', hue_order=hue_order)
plt.savefig('%s/g_value_vs_total_area_ecDNA_sqrt_%s.pdf' % (output_dir, sample))
plt.show()