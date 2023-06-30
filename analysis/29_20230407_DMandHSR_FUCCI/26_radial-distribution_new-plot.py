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

hue_order_r = ['G1', 'S', 'G2']
hue_order = ['unsync', 'G1', 'S', 'G2']
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
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


# new approach
df_r = pd.read_csv('%s%s_radial_cellcycle.txt' % (data_dir, sample), na_values=['.'], sep='\t')

for k in range(len(hue_order)):
    if hue_order[k] == 'unsync':
        data = df
        data_r = df_r
    else:
        data = df[df['cellcycle'] == hue_order[k]].copy().reset_index(drop=True)
        data_r = df_r[df_r['sample'] == hue_order[k]].copy().reset_index(drop=True)
    sns.set_palette(sns.color_palette(line_colors))

    data_r_sort = data_r[data_r['seg'] == 1].copy().reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(9, 6))
    fig.subplots_adjust(left=0.2)
    ax = sns.histplot(data=data_r_sort, x='relative_r', hue='sample_category', hue_order=hue_order1, multiple='dodge', bins=20, weights=data_r_sort['weights'])
    plt.savefig('%s/%s_histplot_DNAFISH_seg_n%s_cell%s_%s.pdf' % (output_dir, figure_name, n_dilation, len(data), hue_order[k]))
    plt.show()

# heatmap
x = np.arange(0.025, 1, 0.05)
x_label = 'relative r'
up = 20
df_seg_normalized = pd.DataFrame(columns=x[:up])
df_int_normalized = pd.DataFrame(columns=x[:up])
df_int_seg_normalized = pd.DataFrame(columns=x[:up])

# seg
for k in range(len(hue_order)):
    if hue_order[k] == 'unsync':
        data = df
        data_r = df_r
    else:
        data = df[df['cellcycle'] == hue_order[k]].copy().reset_index(drop=True)
        data_r = df_r[df_r['sample'] == hue_order[k]].copy().reset_index(drop=True)

    data_heatmap = pd.DataFrame(columns=x[:up])
    for s in hue_order1:
        data_sample = data_r[(data_r['sample_category'] == s) & (data_r['seg'] == 1)].copy().reset_index(drop=True)
        total_data_sample = len(data_sample)
        data_radial = []
        for i in range(len(x[:up])):
            if x[i] == 0.025:
                n_data_radial = len(data_sample[(data_sample['relative_r'] >= 0) & (data_sample['relative_r'] <= 0.05)])
            else:
                n_data_radial = len(data_sample[(data_sample['relative_r'] > (x[i]-0.025)) & (data_sample['relative_r'] <= (x[i]+0.025))])
            data_radial.append(n_data_radial * 1.0/total_data_sample)
        data_heatmap.loc[len(data_heatmap.index)] = data_radial
    data_heatmap.index = hue_order1
    data_heatmap.columns = ['%.3f' % elem for elem in x[:up]]

    plt.subplots(figsize=(12, len(hue_order1)))
    ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
    plt.savefig('%s/%s_heatmap_DNAFISH_seg_before_normalization_n%s_cell%s_%s.pdf' % (output_dir, figure_name, n_dilation, len(data), hue_order[k]))
    plt.show()

    data_heatmap_normalized = pd.DataFrame(columns=x[:up])
    temp = []
    for i in data_heatmap.columns:
        temp.append(data_heatmap[i][hue_order1[1]] / data_heatmap[i][hue_order1[0]])
    data_heatmap_normalized.loc[0] = temp
    data_heatmap_normalized.index = [hue_order1[1]]
    data_heatmap_normalized.columns = ['%.3f' % elem for elem in x[:up]]
    df_seg_normalized.loc[len(df_seg_normalized.index)] = temp

    plt.subplots(figsize=(12, len(hue_order)))
    ax1 = sns.heatmap(data_heatmap_normalized, cbar=0, linewidths=2, vmax=data_heatmap_normalized.values.max(),
                      vmin=data_heatmap_normalized.values.min(), square=True, cmap='coolwarm')
    plt.savefig('%s/%s_heatmap_DNAFISH_seg_n%s_cell%s_%s.pdf' % (output_dir, figure_name, n_dilation, len(data), hue_order[k]))
    plt.show()

# int
for k in range(len(hue_order)):
    if hue_order[k] == 'unsync':
        data = df
        data_r = df_r
    else:
        data = df[df['cellcycle'] == hue_order[k]].copy().reset_index(drop=True)
        data_r = df_r[df_r['sample'] == hue_order[k]].copy().reset_index(drop=True)
    data_heatmap = pd.DataFrame(columns=x[:up])

    for s in hue_order1:
        data_sample = data_r[data_r['sample_category'] == s].copy().reset_index(drop=True)
        data_radial = []
        for i in range(len(x[:up])):
            if x[i] == 0.025:
                filter = (data_sample['relative_r'] >= 0) & (data_sample['relative_r'] <= 0.05)
            else:
                filter = (data_sample['relative_r'] > (x[i]-0.025)) & (data_sample['relative_r'] <= (x[i]+0.025))
            data_radial.append(sum(data_sample[filter]['intensity'])/sum(data_sample['intensity']))
        data_heatmap.loc[len(data_heatmap.index)] = data_radial
    data_heatmap.index = hue_order1
    data_heatmap.columns = ['%.3f' % elem for elem in x[:up]]

    plt.subplots(figsize=(12, len(hue_order1)))
    ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
    plt.savefig('%s/%s_heatmap_DNAFISH_int_before_normalization_n%s_cell%s_%s.pdf' % (output_dir, figure_name, n_dilation, len(data), hue_order[k]))
    plt.show()

    data_heatmap_normalized = pd.DataFrame(columns=x[:up])
    temp = []
    for i in data_heatmap.columns:
        temp.append(data_heatmap[i][hue_order1[1]] / data_heatmap[i][hue_order1[0]])
    data_heatmap_normalized.loc[0] = temp
    data_heatmap_normalized.index = [hue_order1[1]]
    data_heatmap_normalized.columns = ['%.3f' % elem for elem in x[:up]]
    df_int_normalized.loc[len(df_int_normalized.index)] = temp

    plt.subplots(figsize=(12, len(hue_order)))
    ax1 = sns.heatmap(data_heatmap_normalized, cbar=0, linewidths=2, vmax=data_heatmap_normalized.values.max(),
                      vmin=data_heatmap_normalized.values.min(), square=True, cmap='coolwarm')
    plt.savefig('%s/%s_heatmap_DNAFISH_int_n%s_cell%s_%s.pdf' % (output_dir, figure_name, n_dilation, len(data), hue_order[k]))
    plt.show()

# int_seg
for k in range(len(hue_order)):
    if hue_order[k] == 'unsync':
        data = df
        data_r = df_r
    else:
        data = df[df['cellcycle'] == hue_order[k]].copy().reset_index(drop=True)
        data_r = df_r[df_r['sample'] == hue_order[k]].copy().reset_index(drop=True)
    data_heatmap = pd.DataFrame(columns=x[:up])

    for s in hue_order1:
        data_sample = data_r[data_r['sample_category'] == s].copy().reset_index(drop=True)
        data_radial = []
        for i in range(len(x[:up])):
            if x[i] == 0.025:
                filter = (data_sample['relative_r'] >= 0) & (data_sample['relative_r'] <= 0.05)
            else:
                filter = (data_sample['relative_r'] > (x[i]-0.025)) & (data_sample['relative_r'] <= (x[i]+0.025))
            if s == hue_order1[0]:
                feature = 'seg'
            else:
                feature = 'intensity'
            data_radial.append(sum(data_sample[filter][feature])/sum(data_sample[feature]))
        data_heatmap.loc[len(data_heatmap.index)] = data_radial
    data_heatmap.index = hue_order1
    data_heatmap.columns = ['%.3f' % elem for elem in x[:up]]

    plt.subplots(figsize=(12, len(hue_order1)))
    ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
    plt.savefig('%s/%s_heatmap_DNAFISH_int-seg_before_normalization_n%s_cell%s_%s.pdf' % (output_dir, figure_name, n_dilation, len(data), hue_order[k]))
    plt.show()

    data_heatmap_normalized = pd.DataFrame(columns=x[:up])
    temp = []
    for i in data_heatmap.columns:
        temp.append(data_heatmap[i][hue_order1[1]] / data_heatmap[i][hue_order1[0]])
    data_heatmap_normalized.loc[0] = temp
    data_heatmap_normalized.index = [hue_order1[1]]
    data_heatmap_normalized.columns = ['%.3f' % elem for elem in x[:up]]
    df_int_seg_normalized.loc[len(df_int_seg_normalized.index)] = temp

    plt.subplots(figsize=(12, len(hue_order)))
    ax1 = sns.heatmap(data_heatmap_normalized, cbar=0, linewidths=2, vmax=data_heatmap_normalized.values.max(),
                      vmin=data_heatmap_normalized.values.min(), square=True, cmap='coolwarm')
    plt.savefig('%s/%s_heatmap_DNAFISH_int-seg_n%s_cell%s_%s.pdf' % (output_dir, figure_name, n_dilation, len(data), hue_order[k]))
    plt.show()

df_seg_normalized.index = hue_order
df_int_normalized.index = hue_order
df_int_seg_normalized.index = hue_order
df_seg_normalized.columns = ['%.3f' % elem for elem in x[:up]]
df_int_normalized.columns = ['%.3f' % elem for elem in x[:up]]
df_int_seg_normalized.columns = ['%.3f' % elem for elem in x[:up]]

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=df_seg_normalized.values.max(),
                  vmin=df_seg_normalized.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/%s_heatmap_DNAFISH_seg_n%s_cell%s.pdf' % (output_dir, figure_name, n_dilation, len(data)))
plt.show()

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(df_int_normalized, cbar=0, linewidths=2, vmax=df_int_normalized.values.max(),
                  vmin=df_int_normalized.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/%s_heatmap_DNAFISH_int_n%s_cell%s.pdf' % (output_dir, figure_name, n_dilation, len(data)))
plt.show()

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(df_int_seg_normalized, cbar=0, linewidths=2, vmax=df_int_seg_normalized.values.max(),
                  vmin=df_int_seg_normalized.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/%s_heatmap_DNAFISH_int-seg_n%s_cell%s.pdf' % (output_dir, figure_name, n_dilation, len(data)))
plt.show()

plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    if hue_order[k] == 'unsync':
        data = df
    else:
        data = df[df['cellcycle'] == hue_order[k]].copy().reset_index(drop=True)
    plt.plot(x, df_seg_normalized.iloc[k].tolist(), color=line_colors[k], label='%s, n=%s' % (hue_order[k], len(data)))
plt.axhline(y=1, linestyle='--', color='r')
plt.xlabel('relative_r')
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s_rc-new_DNAFISH_seg_n%s.pdf' % (output_dir, figure_name, n_dilation))
plt.show()

plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    if hue_order[k] == 'unsync':
        data = df
    else:
        data = df[df['cellcycle'] == hue_order[k]].copy().reset_index(drop=True)
    plt.plot(x, df_int_normalized.iloc[k].tolist(), color=line_colors[k], label='%s, n=%s' % (hue_order[k], len(data)))
plt.axhline(y=1, linestyle='--', color='r')
plt.xlabel('relative_r')
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s_rc-new_DNAFISH_int_n%s.pdf' % (output_dir, figure_name, n_dilation))
plt.show()

plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    if hue_order[k] == 'unsync':
        data = df
    else:
        data = df[df['cellcycle'] == hue_order[k]].copy().reset_index(drop=True)
    plt.plot(x, df_int_seg_normalized.iloc[k].tolist(), color=line_colors[k], label='%s, n=%s' % (hue_order[k], len(data)))
plt.axhline(y=1, linestyle='--', color='r')
plt.xlabel('relative_r')
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s_rc-new_DNAFISH_int-seg_n%s.pdf' % (output_dir, figure_name, n_dilation))
plt.show()