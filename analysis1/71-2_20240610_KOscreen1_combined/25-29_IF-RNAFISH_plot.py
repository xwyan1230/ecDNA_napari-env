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
import os
import seaborn as sns
from scipy.stats import ttest_ind
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'B2'
reps = [1, 2]
hueorder = ['mCherry', 'GFP']

for rep in reps:
    # 25
    GFP_cutoff = [3.6, 3.0, 3.2, 2.9]
    mCherry_cutoff = [3.1, 3.7, 2.9, 3.4]

    data = pd.read_csv('%s/%s/24_%s_IF_rep%s.txt' % (output_dir, sample, sample, rep), na_values=['.'], sep='\t')

    fig, ax = plt.subplots(figsize=(9, 6))
    fig.subplots_adjust(right=0.8)
    sns.histplot(data=data, x='hoechst')
    plt.show()

    data['log10_GFP'] = np.log10(data['GFP'])
    data['log10_mCherry'] = np.log10(data['mCherry'])
    data['group'] = ['NA'] * len(data)
    data.loc[((data['log10_GFP'] > GFP_cutoff[0]) & (data['log10_mCherry'] < GFP_cutoff[1])) | ((data['log10_GFP'] > GFP_cutoff[2]) & (data['log10_mCherry'] < GFP_cutoff[3])), 'group'] = 'GFP'
    data.loc[((data['log10_mCherry'] > mCherry_cutoff[0]) & (data['log10_GFP'] < mCherry_cutoff[1])) | ((data['log10_mCherry'] > mCherry_cutoff[2]) & (data['log10_GFP'] < mCherry_cutoff[3])), 'group'] = 'mCherry'

    print(data.head())

    plt.subplots(figsize=(9, 9))
    sns.scatterplot(data=data, x='log10_mCherry', y='log10_GFP', alpha=0.5, s=30)
    sns.scatterplot(data=data[data['group'] == 'mCherry'], x='log10_mCherry', y='log10_GFP', color='r', alpha=0.5, s=30)
    sns.scatterplot(data=data[data['group'] == 'GFP'], x='log10_mCherry', y='log10_GFP', color='g', alpha=0.5, s=30)
    plt.savefig('%s/%s/25_%s_IF_rep%s_red-green.pdf' % (output_dir, sample, sample, rep))
    plt.xlim([2.5, 5])
    plt.ylim([3.2, 5])
    plt.savefig('%s/%s/25_%s_IF_rep%s_red-green_fixscale.pdf' % (output_dir, sample, sample, rep))
    plt.show()

    print(len(data[data['group'] == 'GFP']))
    print(len(data[data['group'] == 'mCherry']))

    data.to_csv('%s/%s/25_%s_IF_rep%s_red_green_group.txt' % (output_dir, sample, sample, rep), index=False, sep='\t')

    if os.path.exists("%s/IF_red_green.txt" % output_dir):
        pd_group = pd.read_csv('%s/IF_red_green.txt' % output_dir, na_values=['.'], sep='\t')
    else:
        pd_group = pd.DataFrame(columns=['sample', 'rep', 'GFP_cutoff_0', 'GFP_cutoff_1', 'GFP_cutoff_2', 'GFP_cutoff_3', 'mCherry_cutoff_0', 'mCherry_cutoff_1', 'mCherry_cutoff_2', 'mCherry_cutoff_3'])
    if len(pd_group[(pd_group['sample'] == sample) & (pd_group['rep'] == rep)]) == 0:
        pd_group.loc[len(pd_group.index)] = [sample, rep, GFP_cutoff[0], GFP_cutoff[1], GFP_cutoff[2], GFP_cutoff[3], mCherry_cutoff[0], mCherry_cutoff[1], mCherry_cutoff[2], mCherry_cutoff[3]]
    else:
        location = pd_group[(pd_group['sample'] == sample) & (pd_group['rep'] == rep)].index
        pd_group.loc[location] = [sample, rep, GFP_cutoff[0], GFP_cutoff[1], GFP_cutoff[2], GFP_cutoff[3], mCherry_cutoff[0], mCherry_cutoff[1], mCherry_cutoff[2], mCherry_cutoff[3]]
    pd_group.to_csv('%s/IF_red_green.txt' % output_dir, index=False, sep='\t')

    ### 26
    data = pd.read_csv('%s/%s/25_%s_IF_rep%s_red_green_group.txt' % (output_dir, sample, sample, rep), na_values=['.'],
                       sep='\t')
    data['log10_IF'] = np.log10(data['IF'])
    data_flt = data.sort_values(by='group', ascending=False).copy().reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(4, 6))
    fig.subplots_adjust(left=0.4)
    ax = sns.boxplot(data=data_flt, x='group', y='log10_IF', hue_order=hueorder)
    lines = ax.get_lines()
    categories = ax.get_xticks()
    plt.close()

    fig, ax = plt.subplots(figsize=(4, 6))
    fig.subplots_adjust(left=0.4)
    sns.swarmplot(data=data_flt, x='group', y='log10_IF', hue_order=hueorder, s=1.5,
                  color=(255 / 255, 140 / 255, 0 / 255))

    nobs = [len(data_flt[data_flt['group'] == hueorder[i]]) for i in range(len(hueorder))]
    nobs = ["%s" % str(x) for x in nobs]
    for cat in categories:
        # every 4th line at the interval of 6 is median line
        # 0 -> p25 1 -> p75 2 -> lower whisker 3 -> upper whisker 4 -> p50 5 -> upper extreme value
        y = round(lines[4 + cat * 6].get_ydata()[0], 2)
        ax.text(cat, y, '%s (%s)' % (y, nobs[cat]), ha='center', va='center', fontweight='bold', size=7, color='white',
                bbox=dict(facecolor='#445A64'))

    p = ttest_ind(data_flt[data_flt['group'] == 'mCherry']['log10_IF'].tolist(),
                  data_flt[data_flt['group'] == 'GFP']['log10_IF'].tolist())
    print(p)
    """y, h, col = data_flt['log10_DNAFISH_total_int_merge'].max() + 0.05, 0.05, 'k'
    plt.plot([0, 0, 1, 1], [y, y + h, y + h, y], lw=1, c=col)
    plt.text((0 + 1) * .5, y + h, "ns", ha='center', va='bottom', color=col)"""
    plt.savefig('%s/%s/26_%s_IF_rep%s_MYC.pdf' % (output_dir, sample, sample, rep))
    plt.ylim([2.7, 5])
    plt.savefig('%s/%s/26_%s_IF_rep%s_MYC_fixscale.pdf' % (output_dir, sample, sample, rep))
    plt.show()

    # 28
    GFP_cutoff = [3.6, 3.0, 3.2, 2.85]
    mCherry_cutoff = [3.0, 3.55, 2.85, 3.3]

    data = pd.read_csv('%s/%s/27_%s_RNAFISH_rep%s.txt' % (output_dir, sample, sample, rep), na_values=['.'], sep='\t')
    data['log10_GFP'] = np.log10(data['GFP'])
    data['log10_mCherry'] = np.log10(data['mCherry'])
    data['group'] = ['NA'] * len(data)
    data.loc[((data['log10_GFP'] > GFP_cutoff[0]) & (data['log10_mCherry'] < GFP_cutoff[1])) | (
                (data['log10_GFP'] > GFP_cutoff[2]) & (data['log10_mCherry'] < GFP_cutoff[3])), 'group'] = 'GFP'
    data.loc[((data['log10_mCherry'] > mCherry_cutoff[0]) & (data['log10_GFP'] < mCherry_cutoff[1])) | (
                (data['log10_mCherry'] > mCherry_cutoff[2]) & (
                    data['log10_GFP'] < mCherry_cutoff[3])), 'group'] = 'mCherry'

    print(data.head())

    plt.subplots(figsize=(9, 9))
    sns.scatterplot(data=data, x='log10_mCherry', y='log10_GFP', alpha=0.5, s=30)
    sns.scatterplot(data=data[data['group'] == 'mCherry'], x='log10_mCherry', y='log10_GFP', color='r', alpha=0.5, s=30)
    sns.scatterplot(data=data[data['group'] == 'GFP'], x='log10_mCherry', y='log10_GFP', color='g', alpha=0.5, s=30)
    plt.savefig('%s/%s/28_%s_RNAFISH_rep%s_red-green.pdf' % (output_dir, sample, sample, rep))
    plt.xlim([2.4, 5])
    plt.ylim([3.1, 5])
    plt.savefig('%s/%s/28_%s_RNAFISH_rep%s_red-green_fixscale.pdf' % (output_dir, sample, sample, rep))
    plt.show()

    print(len(data[data['group'] == 'GFP']))
    print(len(data[data['group'] == 'mCherry']))

    data.to_csv('%s/%s/28_%s_RNAFISH_rep%s_red_green_group.txt' % (output_dir, sample, sample, rep), index=False,
                sep='\t')

    if os.path.exists("%s/RNAFISH_red_green.txt" % output_dir):
        pd_group = pd.read_csv('%s/RNAFISH_red_green.txt' % output_dir, na_values=['.'], sep='\t')
    else:
        pd_group = pd.DataFrame(
            columns=['sample', 'rep', 'GFP_cutoff_0', 'GFP_cutoff_1', 'GFP_cutoff_2', 'GFP_cutoff_3',
                     'mCherry_cutoff_0', 'mCherry_cutoff_1', 'mCherry_cutoff_2', 'mCherry_cutoff_3'])
    if len(pd_group[(pd_group['sample'] == sample) & (pd_group['rep'] == rep)]) == 0:
        pd_group.loc[len(pd_group.index)] = [sample, rep, GFP_cutoff[0], GFP_cutoff[1], GFP_cutoff[2], GFP_cutoff[3],
                                             mCherry_cutoff[0], mCherry_cutoff[1], mCherry_cutoff[2], mCherry_cutoff[3]]
    else:
        location = pd_group[(pd_group['sample'] == sample) & (pd_group['rep'] == rep)].index
        pd_group.loc[location] = [sample, rep, GFP_cutoff[0], GFP_cutoff[1], GFP_cutoff[2], GFP_cutoff[3],
                                  mCherry_cutoff[0], mCherry_cutoff[1], mCherry_cutoff[2], mCherry_cutoff[3]]
    pd_group.to_csv('%s/RNAFISH_red_green.txt' % output_dir, index=False, sep='\t')

    ### 29
    data = pd.read_csv('%s/%s/28_%s_RNAFISH_rep%s_red_green_group.txt' % (output_dir, sample, sample, rep),
                       na_values=['.'], sep='\t')
    data['log10_RNAFISH'] = np.log10(data['RNAFISH'])
    data_flt = data.sort_values(by='group', ascending=False).copy().reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(4, 6))
    fig.subplots_adjust(left=0.4)
    ax = sns.boxplot(data=data_flt, x='group', y='log10_RNAFISH', hue_order=hueorder)
    lines = ax.get_lines()
    categories = ax.get_xticks()
    plt.close()

    fig, ax = plt.subplots(figsize=(4, 6))
    fig.subplots_adjust(left=0.4)
    sns.swarmplot(data=data_flt, x='group', y='log10_RNAFISH', hue_order=hueorder, s=1,
                  color=(255 / 255, 140 / 255, 0 / 255))

    nobs = [len(data_flt[data_flt['group'] == hueorder[i]]) for i in range(len(hueorder))]
    nobs = ["%s" % str(x) for x in nobs]
    for cat in categories:
        # every 4th line at the interval of 6 is median line
        # 0 -> p25 1 -> p75 2 -> lower whisker 3 -> upper whisker 4 -> p50 5 -> upper extreme value
        y = round(lines[4 + cat * 6].get_ydata()[0], 2)
        ax.text(cat, y, '%s (%s)' % (y, nobs[cat]), ha='center', va='center', fontweight='bold', size=7, color='white',
                bbox=dict(facecolor='#445A64'))

    p = ttest_ind(data_flt[data_flt['group'] == 'mCherry']['log10_RNAFISH'].tolist(),
                  data_flt[data_flt['group'] == 'GFP']['log10_RNAFISH'].tolist())
    print(p)
    """y, h, col = data_flt['log10_DNAFISH_total_int_merge'].max() + 0.05, 0.05, 'k'
    plt.plot([0, 0, 1, 1], [y, y + h, y + h, y], lw=1, c=col)
    plt.text((0 + 1) * .5, y + h, "ns", ha='center', va='bottom', color=col)"""
    plt.savefig('%s/%s/29_%s_RNAFISH_rep%s_MYC.pdf' % (output_dir, sample, sample, rep))
    plt.ylim([2.7, 5])
    plt.savefig('%s/%s/29_%s_RNAFISH_rep%s_MYC_fixscale.pdf' % (output_dir, sample, sample, rep))
    plt.show()

print("DONE!")