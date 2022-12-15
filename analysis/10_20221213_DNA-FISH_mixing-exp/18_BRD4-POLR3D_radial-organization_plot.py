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

sample = 'DM_mix_DM-H2B-mCherry'
pos_threshold = 10000
neg_threshold = 5000
pos = 'DM H2B-mCherry'
neg = 'DM'
hue_order = [pos, neg]

hist_colors = [(0.90, 0.90, 0.90), (0.95, 0.50, 0.50),  (0.50, 0.90, 0.90)]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25),  (0.30, 0.70, 0.70)]

df = pd.read_csv(("%s%s_ecDNA.txt" % (data_dir, sample)), na_values=['.'], sep='\t')

feature = ['radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_normalized']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

sample_lst = []
for i in range(len(df)):
    if df['mCherry_mean'].tolist()[i] < neg_threshold:
        sample_lst.append(neg)
    elif df['mCherry_mean'].tolist()[i] > pos_threshold:
        sample_lst.append(pos)
    else:
        sample_lst.append('NA')
df['sample'] = sample_lst

# heatmap
data_heatmap = pd.DataFrame(columns=['0.05', '0.15', '0.25', '0.35', '0.45', '0.55', '0.65', '0.75', '0.85', '0.95'])

for s in hue_order:
    data_sample = df[df['sample'] == s].copy().reset_index(drop=True)
    data_radial = pd.DataFrame()
    data_radial['0.05'] = [data_sample['radial_curve_normalized'][i][0] for i in range(len(data_sample))]
    data_radial['0.15'] = [data_sample['radial_curve_normalized'][i][1] for i in range(len(data_sample))]
    data_radial['0.25'] = [data_sample['radial_curve_normalized'][i][2] for i in range(len(data_sample))]
    data_radial['0.35'] = [data_sample['radial_curve_normalized'][i][3] for i in range(len(data_sample))]
    data_radial['0.45'] = [data_sample['radial_curve_normalized'][i][4] for i in range(len(data_sample))]
    data_radial['0.55'] = [data_sample['radial_curve_normalized'][i][5] for i in range(len(data_sample))]
    data_radial['0.65'] = [data_sample['radial_curve_normalized'][i][6] for i in range(len(data_sample))]
    data_radial['0.75'] = [data_sample['radial_curve_normalized'][i][7] for i in range(len(data_sample))]
    data_radial['0.85'] = [data_sample['radial_curve_normalized'][i][8] for i in range(len(data_sample))]
    data_radial['0.95'] = [data_sample['radial_curve_normalized'][i][9] for i in range(len(data_sample))]

    data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()

data_heatmap.index = hue_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_%s.pdf' % (output_dir, sample))
plt.show()

# radial curve
print("Plotting radial curve...")
x = np.arange(0.05, 1, 0.1)
x_label = 'relative r'

for s in hue_order:
    data = df[df['sample'] == s].copy().reset_index(drop=True)
    number_nuclear = len(data)

    mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data['radial_curve_DNAFISH'].tolist())
    mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(data['radial_curve_nuclear'].tolist())
    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_normalized'].tolist())

    plt.subplots(figsize=(12, 9))
    for i in range(len(data)):
        plt.plot(x, data['radial_curve_DNAFISH'][i], alpha=0.05, color=[line_colors[1][j] + 0.05 for j in range(len(line_colors[1]))])
    for i in range(len(data)):
        plt.plot(x, data['radial_curve_nuclear'][i], alpha=0.05, color=[line_colors[0][j]+0.05 for j in range(len(line_colors[0]))])
    plt.plot(x, mean_curve1, color=line_colors[1], label='%s, n=%s' % ('ecDNA (MYC DNA FISH)', number_nuclear))
    plt.plot(x, mean_curve2, color=line_colors[0], label='%s, n=%s' % ('DNA (hoechst stain)', number_nuclear))
    plt.plot(x, ci_lower1, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher1, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_lower2, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher2, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.xlabel(x_label)
    plt.ylim([0.5, 1.5])
    plt.ylabel('radial_curve')
    plt.legend()
    plt.savefig('%sradial_curve_%s_%s.pdf' % (output_dir, sample, s))
    plt.show()

    plt.subplots(figsize=(12, 9))
    for i in range(len(data)):
        plt.plot(x, data['radial_curve_normalized'][i], alpha=0.05, color=[line_colors[1][j] + 0.05 for j in range(len(line_colors[1]))])
    plt.plot(x, mean_curve3, color=line_colors[1], label='%s, n=%s' % ('ecDNA (normalized)', number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.axhline(y=1, color='black', linestyle='--')
    plt.xlabel(x_label)
    plt.ylim([0.5, 1.5])
    plt.ylabel('radial_curve')
    plt.legend()
    plt.savefig('%sradial_curve_normalized_%s_%s.pdf' % (output_dir, sample, s))
    plt.show()

plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['sample'] == hue_order[k]].copy().reset_index(drop=True)
    number_nuclear = len(data)

    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_normalized'].tolist())

    for i in range(len(data)):
        plt.plot(x, data['radial_curve_normalized'][i], alpha=0.05, color=[line_colors[k][j] + 0.05 for j in range(len(line_colors[k]))])
    plt.plot(x, mean_curve3, color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[k], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.5, 1.5])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%sradial_curve_normalized_%s.pdf' % (output_dir, sample))
plt.show()