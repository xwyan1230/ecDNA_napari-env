import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import imutils
import pandas as pd
import shared.image as ima
import tifffile as tif
import math
import nd2
import shared.segmentation as seg
import shared.objects as obj
import shared.dataframe as dat
import os
import seaborn as sns
from scipy.stats import ttest_ind
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'G11'
version = 2
hueorder = ['mCherry', 'GFP']
hue_order1 = ['bg', 'DNAFISH']
total_fov = 16

line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red

### 19
print("Running 19...")
GFP_cutoff = [3.4, 2.9]
mCherry_cutoff = [3, 3.55]

data = pd.read_csv('%s/%s/14_%s_red_green.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
data['log10_GFP'] = np.log10(data['GFP'])
data['log10_mCherry'] = np.log10(data['mCherry'])
data['group'] = ['NA'] * len(data)
data.loc[(data['log10_GFP'] > GFP_cutoff[0]) & (data['log10_mCherry'] < GFP_cutoff[1]), 'group'] = 'GFP'
data.loc[(data['log10_mCherry'] > mCherry_cutoff[0]) & (data['log10_GFP'] < mCherry_cutoff[1]), 'group'] = 'mCherry'

print(data.head())

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data, x='log10_mCherry', y='log10_GFP', alpha=0.5, s=30)
sns.scatterplot(data=data[data['group'] == 'mCherry'], x='log10_mCherry', y='log10_GFP', color='r', alpha=0.5, s=30)
sns.scatterplot(data=data[data['group'] == 'GFP'], x='log10_mCherry', y='log10_GFP', color='g', alpha=0.5, s=30)
plt.savefig('%s/%s/19_%s_red-green.pdf' % (output_dir, sample, sample))
plt.xlim([2.5, 5])
plt.ylim([3.3, 5])
plt.savefig('%s/%s/19_%s_red-green_fixscale.pdf' % (output_dir, sample, sample))
plt.show()

data.to_csv('%s/%s/19_%s_red_green_group.txt' % (output_dir, sample, sample), index=False, sep='\t')

if os.path.exists("%s/red_green.txt" % output_dir):
    pd_group = pd.read_csv('%s/red_green.txt' % output_dir, na_values=['.'], sep='\t')
else:
    pd_group = pd.DataFrame(columns=['sample', 'GFP_cutoff_0', 'GFP_cutoff_1', 'mCherry_cutoff_0', 'mCherry_cutoff_1'])
if len(pd_group[pd_group['sample'] == sample]) == 0:
    pd_group.loc[len(pd_group.index)] = [sample, GFP_cutoff[0], GFP_cutoff[1], mCherry_cutoff[0], mCherry_cutoff[1]]
else:
    location = pd_group[pd_group['sample'] == sample].index
    pd_group.loc[location] = [sample, GFP_cutoff[0], GFP_cutoff[1], mCherry_cutoff[0], mCherry_cutoff[1]]
pd_group.to_csv('%s/red_green.txt' % output_dir, index=False, sep='\t')

### 20
print("Running 20...")
data_rg = pd.read_csv('%s/%s/19_%s_red_green_group.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
data = pd.read_csv('%s/%s/16_%s_copy_number.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
data = data.sort_values(['sample', 'fov', 'label_mean_int'], ascending=[True, True, True]).copy().reset_index(drop=True)
data_rg = data_rg.sort_values(['sample', 'fov', 'label_mean_int'], ascending=[True, True, True]).copy().reset_index(drop=True)

if (data['label_mean_int'].tolist() == data_rg['label_mean_int'].tolist()) & (data['fov'].tolist() == data_rg['fov'].tolist()):
    data['group'] = data_rg['group']
    data_flt = data[data['group'].isin(['GFP', 'mCherry'])].copy().reset_index(drop=True)
    data_flt['log10_DNAFISH_total_int_merge'] = np.log10(data_flt['DNAFISH_total_int_merge']+1)
    data_flt = data_flt.sort_values(by='group', ascending=False).copy().reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(4, 6))
    fig.subplots_adjust(left=0.4)
    ax = sns.boxplot(data=data_flt, x='group', y='log10_DNAFISH_total_int_merge', hue_order=hueorder)
    lines = ax.get_lines()
    categories = ax.get_xticks()
    plt.close()

    fig, ax = plt.subplots(figsize=(4, 6))
    fig.subplots_adjust(left=0.4)
    sns.swarmplot(data=data_flt, x='group', y='log10_DNAFISH_total_int_merge', hue_order=hueorder, s=2, color=(255 / 255, 140 / 255, 0 / 255))

    nobs = [len(data_flt[data_flt['group'] == hueorder[i]]) for i in range(len(hueorder))]
    nobs = ["%s" % str(x) for x in nobs]
    for cat in categories:
        # every 4th line at the interval of 6 is median line
        # 0 -> p25 1 -> p75 2 -> lower whisker 3 -> upper whisker 4 -> p50 5 -> upper extreme value
        y = round(lines[4 + cat * 6].get_ydata()[0], 2)
        ax.text(cat, y, '%s (%s)' % (y, nobs[cat]), ha='center', va='center', fontweight='bold', size=7, color='white', bbox=dict(facecolor='#445A64'))

    p = ttest_ind(data_flt[data_flt['group'] == 'mCherry']['DNAFISH_total_int_merge'].tolist(), data_flt[data_flt['group'] == 'GFP']['DNAFISH_total_int_merge'].tolist())
    print(p)
    """y, h, col = data_flt['log10_DNAFISH_total_int_merge'].max() + 0.05, 0.05, 'k'
    plt.plot([0, 0, 1, 1], [y, y + h, y + h, y], lw=1, c=col)
    plt.text((0 + 1) * .5, y + h, "ns", ha='center', va='bottom', color=col)"""
    plt.savefig('%s/%s/20_%s_copy_number.pdf' % (output_dir, sample, sample))
    plt.ylim([7, 9])
    plt.savefig('%s/%s/20_%s_copy_number_fixscale.pdf' % (output_dir, sample, sample))
    plt.show()

    data.to_csv('%s/%s/20_%s_copy_number_group.txt' % (output_dir, sample, sample), index=False, sep='\t')

else:
    print("no!")

### 21
print("Running 21...")
data_rg = pd.read_csv('%s/%s/19_%s_red_green_group.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
if version == 1:
    data = pd.read_csv('%s/%s/17_%s_cluster.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
else:
    data = pd.read_csv('%s/%s/31_%s_cluster.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')

data = data.sort_values(['sample', 'fov', 'label_nuclear'], ascending=[True, True, True]).copy().reset_index(drop=True)
data_rg = data_rg.sort_values(['sample', 'fov', 'label_mean_int'], ascending=[True, True, True]).copy().reset_index(
    drop=True)

if (data['label_nuclear'].tolist() == data_rg['label_mean_int'].tolist()) & (
        data['fov'].tolist() == data_rg['fov'].tolist()):
    data['group'] = data_rg['group']
    data_flt = data[(data['group'].isin(['GFP', 'mCherry'])) & (data['n_ecDNA'] != -1) & (data['total_area_ratio_ecDNA'] < 0.5)].copy().reset_index(drop=True)
    data_flt = data_flt.sort_values(by='group', ascending=False).copy().reset_index(drop=True)
    feature = ['area_ind_ecDNA', 'area_ratio_ind_ecDNA', 'percentage_area_curve_ecDNA', 'cum_area_ind_ecDNA_filled']
    for f in feature:
        data_flt[f] = [dat.str_to_float(data_flt[f][i]) for i in range(len(data_flt))]
    data_flt1 = data_flt[data_flt['total_area_ratio_ecDNA'] > 0.02].copy().reset_index(drop=True)

    features = ['n_ecDNA', 'total_area_ecDNA', 'total_area_ratio_ecDNA']
    f_range = [[-0.5, 25], [-50, 2000], [-0.01, 0.3]]

    for i in range(len(features)):
        feature = features[i]
        fig, ax = plt.subplots(figsize=(4, 6))
        fig.subplots_adjust(left=0.4)
        ax = sns.boxplot(data=data_flt, x='group', y=feature, hue_order=hueorder)
        lines = ax.get_lines()
        categories = ax.get_xticks()
        plt.close()

        fig, ax = plt.subplots(figsize=(4, 6))
        fig.subplots_adjust(left=0.4)
        sns.swarmplot(data=data_flt, x='group', y=feature, hue_order=hueorder, s=1,
                      color=(255 / 255, 140 / 255, 0 / 255))

        nobs = [len(data_flt[data_flt['group'] == hueorder[i]]) for i in range(len(hueorder))]
        nobs = ["%s" % str(x) for x in nobs]
        for cat in categories:
            # every 4th line at the interval of 6 is median line
            # 0 -> p25 1 -> p75 2 -> lower whisker 3 -> upper whisker 4 -> p50 5 -> upper extreme value
            y = round(lines[4 + cat * 6].get_ydata()[0], 2)
            ax.text(cat, y, '%s (%s)' % (y, nobs[cat]), ha='center', va='center', fontweight='bold', size=7,
                    color='white', bbox=dict(facecolor='#445A64'))

        p = ttest_ind(data_flt[data_flt['group'] == 'mCherry'][feature].tolist(),
                      data_flt[data_flt['group'] == 'GFP'][feature].tolist())
        print(p)
        """y, h, col = data_flt['log10_DNAFISH_total_int_merge'].max() + 0.05, 0.05, 'k'
        plt.plot([0, 0, 1, 1], [y, y + h, y + h, y], lw=1, c=col)
        plt.text((0 + 1) * .5, y + h, "ns", ha='center', va='bottom', color=col)"""
        if version == 1:
            plt.savefig('%s/%s/21_%s_%s.pdf' % (output_dir, sample, sample, feature))
            plt.ylim(f_range[i])
            plt.savefig('%s/%s/21_%s_%s_fixscale.pdf' % (output_dir, sample, sample, feature))
            plt.show()
        else:
            plt.savefig('%s/%s/33_%s_%s.pdf' % (output_dir, sample, sample, feature))
            plt.ylim(f_range[i])
            plt.savefig('%s/%s/33_%s_%s_fixscale.pdf' % (output_dir, sample, sample, feature))
            plt.show()

    features = ['percentage_area_curve_ecDNA', 'cum_area_ind_ecDNA_filled']
    for f in features:
        x_label = 'number of ecDNA hub'
        df1 = data_flt1[data_flt1['group'] == hueorder[0]].copy().reset_index(drop=True)
        df2 = data_flt1[data_flt1['group'] == hueorder[1]].copy().reset_index(drop=True)
        number_nuclear1 = len(df1)
        number_nuclear2 = len(df2)
        mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(dat.list_fill_with_last_num(df1[f].tolist()))
        mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(dat.list_fill_with_last_num(df2[f].tolist()))

        plt.subplots(figsize=(6, 4))
        for i in range(len(df1)):
            plt.plot(df1[f][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
        for i in range(len(df2)):
            plt.plot(df2[f][i], alpha=0.05, color=[line_colors[1][j] + 0.05 for j in range(len(line_colors[1]))])
        plt.plot(mean_curve1, color=line_colors[0], label='%s, n=%s' % (hueorder[0], number_nuclear1))
        plt.plot(mean_curve2, color=line_colors[1], label='%s, n=%s' % (hueorder[1], number_nuclear2))
        plt.plot(ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
        plt.plot(ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
        plt.plot(ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
        plt.plot(ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
        plt.xlabel(x_label)
        plt.ylabel(f)
        plt.legend()
        if version == 1:
            plt.savefig('%s/%s/21_%s_%s.pdf' % (output_dir, sample, sample, f))
            plt.xlim([-0.5, 25])
            if f == 'cum_area_ind_ecDNA_filled':
                plt.ylim([-50, 2000])
            plt.savefig('%s/%s/21_%s_%s_fixscale.pdf' % (output_dir, sample, sample, f))
            plt.show()
        else:
            plt.savefig('%s/%s/33_%s_%s.pdf' % (output_dir, sample, sample, f))
            plt.xlim([-0.5, 25])
            if f == 'cum_area_ind_ecDNA_filled':
                plt.ylim([-50, 2000])
            plt.savefig('%s/%s/33_%s_%s_fixscale.pdf' % (output_dir, sample, sample, f))
            plt.show()

    features = ['per_AUC']
    f_range = [[0, 1]]

    for i in range(len(features)):
        feature = features[i]
        fig, ax = plt.subplots(figsize=(4, 6))
        fig.subplots_adjust(left=0.4)
        ax = sns.boxplot(data=data_flt1, x='group', y=feature, hue_order=hueorder)
        lines = ax.get_lines()
        categories = ax.get_xticks()
        plt.close()

        fig, ax = plt.subplots(figsize=(4, 6))
        fig.subplots_adjust(left=0.4)
        sns.swarmplot(data=data_flt1, x='group', y=feature, hue_order=hueorder, s=1,
                      color=(255 / 255, 140 / 255, 0 / 255))

        nobs = [len(data_flt1[data_flt1['group'] == hueorder[i]]) for i in range(len(hueorder))]
        nobs = ["%s" % str(x) for x in nobs]
        for cat in categories:
            # every 4th line at the interval of 6 is median line
            # 0 -> p25 1 -> p75 2 -> lower whisker 3 -> upper whisker 4 -> p50 5 -> upper extreme value
            y = round(lines[4 + cat * 6].get_ydata()[0], 2)
            ax.text(cat, y, '%s (%s)' % (y, nobs[cat]), ha='center', va='center', fontweight='bold', size=7,
                    color='white', bbox=dict(facecolor='#445A64'))

        p = ttest_ind(data_flt1[data_flt1['group'] == 'mCherry'][feature].tolist(),
                      data_flt1[data_flt1['group'] == 'GFP'][feature].tolist())
        print(p)
        """y, h, col = data_flt['log10_DNAFISH_total_int_merge'].max() + 0.05, 0.05, 'k'
        plt.plot([0, 0, 1, 1], [y, y + h, y + h, y], lw=1, c=col)
        plt.text((0 + 1) * .5, y + h, "ns", ha='center', va='bottom', color=col)"""
        if version == 1:
            plt.savefig('%s/%s/21_%s_%s.pdf' % (output_dir, sample, sample, feature))
            plt.ylim(f_range[i])
            plt.savefig('%s/%s/21_%s_%s_fixscale.pdf' % (output_dir, sample, sample, feature))
            plt.show()
        else:
            plt.savefig('%s/%s/33_%s_%s.pdf' % (output_dir, sample, sample, feature))
            plt.ylim(f_range[i])
            plt.savefig('%s/%s/33_%s_%s_fixscale.pdf' % (output_dir, sample, sample, feature))
            plt.show()
    if version == 1:
        data.to_csv('%s/%s/21_%s_cluster_group.txt' % (output_dir, sample, sample), index=False, sep='\t')
    else:
        data.to_csv('%s/%s/33_%s_cluster_group.txt' % (output_dir, sample, sample), index=False, sep='\t')

else:
    print('no')

### 22
print("Running 22...")
data_rg = pd.read_csv('%s/%s/19_%s_red_green_group.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
if version == 1:
    data = pd.read_csv('%s/%s/18_%s_radial.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
else:
    data = pd.read_csv('%s/%s/32_%s_radial.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')

data = data.sort_values(['sample', 'fov', 'label_nuclear'], ascending=[True, True, True]).copy().reset_index(drop=True)
data_rg = data_rg.sort_values(['sample', 'fov', 'label_mean_int'], ascending=[True, True, True]).copy().reset_index(drop=True)

if version == 1:
    pd_thresh = pd.read_csv('%s/%s/17_%s_DNAFISH_thresh.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
else:
    pd_thresh = pd.read_csv('%s/%s/31_%s_DNAFISH_thresh.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
fov_flt = pd_thresh[pd_thresh['DNAFISH_thresh'] > 1000]['fov'].tolist()

if (data['label_nuclear'].tolist() == data_rg['label_mean_int'].tolist()) & (data['fov'].tolist() == data_rg['fov'].tolist()):
    data['group'] = data_rg['group']
    data_flt = data[(data['group'].isin(['GFP', 'mCherry'])) & (data['fov'].isin(fov_flt))].copy().reset_index(drop=True)
    feature_lst = ['DNAFISH_seg_label', 'int_r_to_edge', 'int_relative_r']
    for f in feature_lst:
        data_flt[f] = [dat.str_to_float(data_flt[f][i]) for i in range(len(data_flt))]

    df_r = pd.DataFrame()
    for h in range(len(hueorder)):
        df = data_flt[data_flt['group'] == hueorder[h]].copy().reset_index(drop=True)
        # new approach
        df_r_temp = pd.DataFrame()
        seg_lst = []
        relative_r_lst = []
        absolute_r_lst = []
        for i in range(len(df)):
            relative_r_lst = relative_r_lst + df['int_relative_r'][i]
            seg_lst = seg_lst + df['DNAFISH_seg_label'][i]
            absolute_r_lst = absolute_r_lst + df['int_r_to_edge'][i]
        df_r_temp['relative_r'] = relative_r_lst
        df_r_temp['seg'] = seg_lst
        df_r_temp['absolute_r'] = absolute_r_lst
        df_r_temp['group'] = [hueorder[h]] * len(df_r_temp)
        len_sample = len(df_r_temp[df_r_temp['seg'] == 1.0])
        print(len_sample)
        # df_r_temp['weights'] = [1.0 / len_sample] * len(df_r_temp)
        df_r = pd.concat([df_r, df_r_temp], axis=0)
    if version == 1:
        df_r.to_csv('%s/%s/22_%s_radial_calculated.txt' % (output_dir, sample, sample), index=False, sep='\t')
    else:
        df_r.to_csv('%s/%s/34_%s_radial_calculated.txt' % (output_dir, sample, sample), index=False, sep='\t')

### 23
print("Running 23...")
if version == 1:
    data = pd.read_csv('%s/%s/22_%s_radial_calculated.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
else:
    data = pd.read_csv('%s/%s/34_%s_radial_calculated.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
df = data[data['group'].isin(['GFP', 'mCherry'])].copy().reset_index(drop=True)

# heatmap
x = np.arange(0.0125, 1, 0.025)
x_label = 'relative r'
up = 40
df_seg_normalized = pd.DataFrame(columns=x[:up])

# seg
for k in range(len(hueorder)):
    data_r = df[(df['group'] == hueorder[k])].copy().reset_index(drop=True)
    data_r_sort = data_r[data_r['seg'] == 1].copy().reset_index(drop=True)
    data_r_sort['sample_category'] = ['DNAFISH'] * len(data_r_sort)
    data_r_bg = data_r.copy()
    data_r_bg['seg'] = [1] * len(data_r_bg)
    data_r_bg['sample_category'] = ['bg'] * len(data_r_bg)
    data_r_final = pd.concat([data_r_sort, data_r_bg], axis=0)

    data_heatmap = pd.DataFrame(columns=x[:up])
    for s in hue_order1:
        data_sample = data_r_final[(data_r_final['sample_category'] == s) & (data_r_final['seg'] == 1)].copy().reset_index(drop=True)
        total_data_sample = len(data_sample)
        data_radial = []
        for i in range(len(x[:up])):
            if x[i] == 0.0125:
                n_data_radial = len(data_sample[(data_sample['relative_r'] >= 0) & (data_sample['relative_r'] <= 0.025)])
            else:
                n_data_radial = len(data_sample[(data_sample['relative_r'] > (x[i]-0.0125)) & (data_sample['relative_r'] <= (x[i]+0.0125))])
            data_radial.append(n_data_radial * 1.0/total_data_sample)
        data_heatmap.loc[len(data_heatmap.index)] = data_radial
    data_heatmap.index = hue_order1
    data_heatmap.columns = ['%.3f' % elem for elem in x[:up]]

    temp = []
    for i in data_heatmap.columns:
        temp.append(data_heatmap[i][hue_order1[1]] / data_heatmap[i][hue_order1[0]])
    df_seg_normalized.loc[len(df_seg_normalized.index)] = temp
df_seg_normalized.index = hueorder
df_seg_normalized.columns = ['%.3f' % elem for elem in x[:up]]

if version == 1:
    df_seg_normalized.to_csv('%s%s/23_%s_radial_summary_relativer.txt' % (output_dir, sample, sample), index=False, sep='\t')
else:
    df_seg_normalized.to_csv('%s%s/35_%s_radial_summary_relativer.txt' % (output_dir, sample, sample), index=False,
                             sep='\t')

plt.subplots(figsize=(12, len(hueorder)))
ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=df_seg_normalized.values.max(),
                  vmin=df_seg_normalized.values.min(), square=True, cmap='coolwarm')
if version == 1:
    plt.savefig('%s%s/23_%s_heatmap_DNAFISH_seg.pdf' % (output_dir, sample, sample))
else:
    plt.savefig('%s%s/35_%s_heatmap_DNAFISH_seg.pdf' % (output_dir, sample, sample))
plt.close()

plt.subplots(figsize=(12, len(hueorder)))
ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=1.8,
                  vmin=0.2, square=True, cmap='coolwarm')
if version == 1:
    plt.savefig('%s%s/23_%s_heatmap_DNAFISH_seg_fixscale.pdf' % (output_dir, sample, sample))
else:
    plt.savefig('%s%s/35_%s_heatmap_DNAFISH_seg_fixscale.pdf' % (output_dir, sample, sample))
plt.show()

# radial_curve
x = np.arange(0.0125, 1, 0.025)
x_label = 'relative r'
plt.subplots(figsize=(12, 9))
for k in range(len(hueorder)):
    mean_curve = df_seg_normalized.iloc[k]
    plt.plot(x, mean_curve, color=line_colors[k], label='%s' % hueorder[k])
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
# plt.ylim([0.4, 1.6])
plt.ylabel('radial_curve')
plt.legend()
if version == 1:
    plt.savefig('%s%s/23_%s_radial_curve_DNAFISH_seg.pdf' % (output_dir, sample, sample))
    plt.ylim([0.2, 1.8])
    plt.savefig('%s%s/23_%s_radial_curve_DNAFISH_seg_fixscale.pdf' % (output_dir, sample, sample))
    plt.show()
else:
    plt.savefig('%s%s/35_%s_radial_curve_DNAFISH_seg.pdf' % (output_dir, sample, sample))
    plt.ylim([0, 2])
    plt.savefig('%s%s/35_%s_radial_curve_DNAFISH_seg_fixscale.pdf' % (output_dir, sample, sample))
    plt.show()

# heatmap
# x = np.arange(98.25, 0, -2.5)
x = np.arange(39.5, 0, -1)
x_label = 'absolute r'
up = 40
df_seg_normalized = pd.DataFrame(columns=x[:up])

# seg
for k in range(len(hueorder)):
    data_r = df[(df['group'] == hueorder[k])].copy().reset_index(drop=True)
    data_r_sort = data_r[data_r['seg'] == 1].copy().reset_index(drop=True)
    data_r_sort['sample_category'] = ['DNAFISH'] * len(data_r_sort)
    data_r_bg = data_r.copy()
    data_r_bg['seg'] = [1] * len(data_r_bg)
    data_r_bg['sample_category'] = ['bg'] * len(data_r_bg)
    data_r_final = pd.concat([data_r_sort, data_r_bg], axis=0)

    data_heatmap = pd.DataFrame(columns=x[:up])
    for s in hue_order1:
        data_sample = data_r_final[(data_r_final['sample_category'] == s) & (data_r_final['seg'] == 1)].copy().reset_index(drop=True)
        total_data_sample = len(data_sample)
        data_radial = []
        for i in range(len(x[:up])):
            if x[i] == 0.75:
                n_data_radial = len(data_sample[(data_sample['absolute_r'] >= 0) & (data_sample['absolute_r'] <= 1.5)])
            else:
                n_data_radial = len(data_sample[(data_sample['absolute_r'] > (x[i]-0.75)) & (data_sample['absolute_r'] <= (x[i]+0.75))])
            data_radial.append(n_data_radial * 1.0/total_data_sample)
        data_heatmap.loc[len(data_heatmap.index)] = data_radial
    data_heatmap.index = hue_order1
    data_heatmap.columns = ['%.3f' % elem for elem in x[:up]]

    temp = []
    for i in data_heatmap.columns:
        temp.append(data_heatmap[i][hue_order1[1]] / data_heatmap[i][hue_order1[0]])
    df_seg_normalized.loc[len(df_seg_normalized.index)] = temp
df_seg_normalized.index = hueorder
df_seg_normalized.columns = ['%.3f' % elem for elem in x[:up]]

if version == 1:
    df_seg_normalized.to_csv('%s%s/23_%s_radial_summary_absoluter.txt' % (output_dir, sample, sample), index=False, sep='\t')
else:
    df_seg_normalized.to_csv('%s%s/35_%s_radial_summary_absoluter.txt' % (output_dir, sample, sample), index=False,
                             sep='\t')

plt.subplots(figsize=(12, len(hueorder)))
ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=df_seg_normalized.values.max(),
                  vmin=df_seg_normalized.values.min(), square=True, cmap='coolwarm')
if version == 1:
    plt.savefig('%s%s/23_%s_heatmap_edge_DNAFISH_seg.pdf' % (output_dir, sample, sample))
else:
    plt.savefig('%s%s/35_%s_heatmap_edge_DNAFISH_seg.pdf' % (output_dir, sample, sample))
plt.close()

plt.subplots(figsize=(12, len(hueorder)))
ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=1.8,
                  vmin=0.2, square=True, cmap='coolwarm')
if version == 1:
    plt.savefig('%s%s/23_%s_heatmap_edge_DNAFISH_seg_fixscale.pdf' % (output_dir, sample, sample))
else:
    plt.savefig('%s%s/35_%s_heatmap_edge_DNAFISH_seg_fixscale.pdf' % (output_dir, sample, sample))
plt.show()

# radial_curve
x = np.arange(1.25, 100, 2.5)
x_label = 'absolute r'
plt.subplots(figsize=(12, 9))
for k in range(len(hueorder)):
    mean_curve = df_seg_normalized.iloc[k]
    plt.plot(x, mean_curve, color=line_colors[k], label='%s' % hueorder[k])
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
# plt.ylim([0.4, 1.6])
plt.ylabel('radial_curve')
plt.legend()
if version == 1:
    plt.savefig('%s%s/23_%s_radial_curve_edge_DNAFISH_seg.pdf' % (output_dir, sample, sample))
    plt.ylim([0.2, 1.8])
    plt.savefig('%s%s/23_%s_radial_curve_edge_DNAFISH_seg_fixscale.pdf' % (output_dir, sample, sample))
    plt.show()
else:
    plt.savefig('%s%s/35_%s_radial_curve_edge_DNAFISH_seg.pdf' % (output_dir, sample, sample))
    plt.ylim([0, 2])
    plt.savefig('%s%s/35_%s_radial_curve_edge_DNAFISH_seg_fixscale.pdf' % (output_dir, sample, sample))
    plt.show()

print("DONE!")