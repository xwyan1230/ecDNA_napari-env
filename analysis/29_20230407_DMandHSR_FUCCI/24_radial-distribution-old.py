import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import numpy as np
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

hue_order = ['unsync', 'G1', 'S', 'G2']
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
line_colors1 = [(119/255, 136/255, 153/255), (220/255, 20/255, 60/255), (154/255, 205/255, 50/255), (255/255, 127/255, 80/255)]
hue_order1 = ['background', 'DNAFISH']

feature = ['radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_normalized', 'radial_curve_DNAFISH_seg',
           'nuclear_int', 'DNAFISH_int', 'DNAFISH_seg_label', 'int_r_to_edge', 'int_r_to_center', 'int_relative_r']

for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

df_sample = df[df['cellcycle'].isin(hue_order)].copy().reset_index(drop=True)
df = df_sample
print(len(df))


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask) if e != 0]
    out = list(map(img.__getitem__, index))
    return out


# heatmap
column_lst = ['0.025', '0.075', '0.125', '0.175', '0.225', '0.275', '0.325', '0.375', '0.425', '0.475', '0.525',
              '0.575', '0.625', '0.675', '0.725', '0.775', '0.825', '0.875', '0.925', '0.975']

data_heatmap = pd.DataFrame(columns=column_lst)

for s in hue_order:
    if s == 'unsync':
        data_sample = df
    else:
        data_sample = df[df['cellcycle'] == s].copy().reset_index(drop=True)
    data_radial = pd.DataFrame()
    for j in range(len(column_lst)):
        data_radial[column_lst[j]] = [data_sample['radial_curve_DNAFISH'][i][j] for i in range(len(data_sample))]

    data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()

data_heatmap.index = hue_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_DNAFISH_n%s_%s.pdf' % (output_dir, n_dilation, figure_name))
plt.show()

# radial curve
print("Plotting radial curve...")
x = np.arange(0.025, 1, 0.05)
x_label = 'relative r'

for s in hue_order:
    if s == 'unsync':
        data = df
    else:
        data = df[df['cellcycle'] == s].copy().reset_index(drop=True)
    number_nuclear = len(data)

    mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data['radial_curve_DNAFISH'].tolist())
    mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(data['radial_curve_nuclear'].tolist())

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
    plt.savefig('%sradial_curve_n%s_%s_%s.pdf' % (output_dir, n_dilation, figure_name, s))
    plt.show()

plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    if hue_order[k] == 'unsync':
        data = df
    else:
        data = df[df['cellcycle'] == hue_order[k]].copy().reset_index(drop=True)
    number_nuclear = len(data)

    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_DNAFISH'].tolist())

    for i in range(len(data)):
        plt.plot(x, data['radial_curve_DNAFISH'][i], alpha=0.01, color=line_colors1[k])
    plt.plot(x, mean_curve3, color=line_colors1[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors1[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors1[k], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.5, 1.5])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%sradial_curve_DNAFISH_n%s_%s.pdf' % (output_dir, n_dilation, figure_name))
plt.show()

