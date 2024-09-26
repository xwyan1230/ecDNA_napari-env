import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import skimage.io as skio
import math
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
data_dir1 = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

hc = [5000, 40000]
cutoff = 3.2

folder = 'point5uM_48hr'
batch = 'point5uM_48hr_3'
if int(batch[-1]) == 1:
    samples = ['XY02', 'XY04', 'XY06', 'XY08', 'XY10', 'XY12', 'XY14', 'XY16', 'XY18', 'XY20']
    samples = ['XY96']
    # samples = ['XY0%s' % (x+1) for x in range(9)] + ['XY%s' % (x+10) for x in range(191)]
elif int(batch[-1]) == 2:
    samples = ['XY%s' % (2*x+111) for x in range(5)] + ['XY%s' % (2*x+101) for x in range(5)]
    samples = ['XY192']
    # samples = ['XY70', 'XY80', 'XY90', 'XY150', 'XY170']
    # samples = ['XY02', 'XY04', 'XY06', 'XY08', 'XY10', 'XY12', 'XY14', 'XY16', 'XY18', 'XY20']
else:
    samples = ['XY01', 'XY02']+['XY%s' % (2*x+201) for x in range(5)] + ['XY%s' % (2*x+202) for x in range(5)]
    samples = ['XY202']
    # samples = ['XY0%s' % (x+1) for x in range(9)] + ['XY%s' % (x+10) for x in range(5)]

identities = ['neg'] * 5 + ['pos'] * 5 + ['mix'] * 40
# colors = ['#808080', '#bc4d4a'] * 5 + []

df = pd.DataFrame(columns=['sample', 'identity', 'G1', 'S', 'G2M'])
for i in range(len(samples)):
    sample = samples[i]
    print(sample)
    # identity = identities[i]
    data = pd.read_csv('%s/%s/%s/txt/%s.txt' % (data_dir, folder, batch, sample), na_values=['.'], sep='\t')
    data['log10_emiRFP670'] = np.log10(data['emiRFP670'])
    data['log2_H2-3'] = np.log2(data['H2-3'])
    data['log2_AzaleaB5'] = np.log2(data['AzaleaB5'])

    plt.subplots(figsize=(9, 7))
    m = plt.hist(data['log2_H2-3'], bins=50)
    red_index = m[0].tolist().index(np.max(m[0][:-2]))
    red_axis = m[1][red_index]
    if red_index >= 15:
        red_axis = 10.1
    print(red_index)
    print(red_axis)
    plt.close()

    plt.subplots(figsize=(9, 7))
    m = plt.hist(data['log2_AzaleaB5'], bins=50)
    green_index = m[0].tolist().index(np.max(m[0][:-2]))
    green_axis = m[1][green_index]
    if green_index >= 15:
        green_axis = 8.3
    print(green_index)
    print(green_axis)
    plt.close()

    green_cutoff = green_axis+0.65
    red_cutoff = red_axis+0.65
    k = (16-red_cutoff)/(16-green_cutoff)
    b = (red_cutoff-green_cutoff)*16/(16-green_cutoff)

    data = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])].reset_index(drop=True)
    data_G1 = data[(data['log2_AzaleaB5'] > green_cutoff) & (data['log2_H2-3'] <= red_cutoff)]
    data_G1S = data[(data['log2_AzaleaB5'] <= green_cutoff) & (data['log2_H2-3'] < red_cutoff)]
    data_S = data[(data['log2_AzaleaB5'] < green_cutoff) & (data['log2_H2-3'] >= red_cutoff)]
    data_G2M = data[(data['log2_AzaleaB5'] >= green_cutoff) & (data['log2_H2-3']-k*data['log2_AzaleaB5']-b > 0)]
    data_G2MG1 = data[(data['log2_H2-3'] > red_cutoff) & (data['log2_H2-3']-k*data['log2_AzaleaB5']-b <= 0)]

    fig, ax = plt.subplots(figsize=(9, 7))
    # fig.subplots_adjust(right=0.8)
    sns.scatterplot(data=data, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#cccccc') ##4d4e4e
    sns.scatterplot(data=data_G1, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#bc4d4a')
    sns.scatterplot(data=data_G1S, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#3d5a80')
    sns.scatterplot(data=data_S, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#669daa')
    sns.scatterplot(data=data_G2M, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#dd933c')
    sns.scatterplot(data=data_G2MG1, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#d2b48c')
    plt.xlim([7, 16.5])
    plt.ylim([9, 16.5])
    # plt.savefig('%s%s/%s_cell-cycle_%s.pdf' % (output_dir, batch, batch, sample))
    plt.show()

    """plt.subplots(figsize=(9, 7))
    plt.hist([data_G1['hoechst'], data_S['hoechst'], data_G2M['hoechst']], weights=[np.ones(len(data_G1)) / len(data_G1), np.ones(len(data_S)) / len(data_S), np.ones(len(data_G2M)) / len(data_G2M)],
             range=[7000, 42000], bins=30, color=['#bc4d4a','#669daa', '#dd933c'], edgecolor='w')
    # plt.savefig('%s/%s/%s_cell-cycle_hoechst_%s.pdf' % (output_dir, batch, batch, sample))
    plt.show()"""

df.to_csv('%s/%s/%s_cell-cycle.txt' % (data_dir, folder, batch), index=False, sep='\t')

