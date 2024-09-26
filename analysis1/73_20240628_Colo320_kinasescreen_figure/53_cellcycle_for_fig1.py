import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import skimage.io as skio
import os
import math

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
data_dir1 = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

hc = [5000, 40000]
cutoff = 3.2

batch = '24hr_density'
samples = ['XY%s' % (x+201) for x in range(10)]
identities = ['neg', 'pos'] * 5 + ['mix'] * 10
samples =['XY%s' % (x+201) for x in range(10)] + ['XY%s' % (x+121) for x in range(10)]


df = pd.DataFrame(columns=['sample', 'identity', 'G1', 'G1S', 'S', 'G2M', 'G2MG1'])
for i in range(len(samples)):
    sample = samples[i]
    identity = identities[i]
    data = pd.read_csv('%s/%s/%s/txt/%s.txt' % (data_dir, batch, batch, sample), na_values=['.'], sep='\t')
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

    green_cutoff = green_axis + 0.65
    red_cutoff = red_axis + 0.65
    k = (16 - red_cutoff) / (16 - green_cutoff)
    b = (red_cutoff - green_cutoff) * 16 / (16 - green_cutoff)

    data = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])].reset_index(drop=True)
    data_G1 = data[(data['log2_AzaleaB5'] > green_cutoff) & (data['log2_H2-3'] <= red_cutoff)]
    data_G1S = data[(data['log2_AzaleaB5'] <= green_cutoff) & (data['log2_H2-3'] < red_cutoff)]
    data_S = data[(data['log2_AzaleaB5'] < green_cutoff) & (data['log2_H2-3'] >= red_cutoff)]
    data_G2M = data[(data['log2_AzaleaB5'] >= green_cutoff) & (data['log2_H2-3'] - k * data['log2_AzaleaB5'] - b > 0)]
    data_G2MG1 = data[(data['log2_H2-3'] > red_cutoff) & (data['log2_H2-3'] - k * data['log2_AzaleaB5'] - b <= 0)]

    fig, ax = plt.subplots(figsize=(9, 7))
    # fig.subplots_adjust(right=0.8)
    sns.scatterplot(data=data, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#cccccc')  ##4d4e4e
    sns.scatterplot(data=data_G1, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#bc4d4a') # '#bc4d4a', '#3d5a80', '#669daa', '#dd933c', '#d2b48c'
    sns.scatterplot(data=data_G1S, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#3d5a80')
    sns.scatterplot(data=data_S, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#669daa')
    sns.scatterplot(data=data_G2M, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#dd933c')
    sns.scatterplot(data=data_G2MG1, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#d2b48c')
    plt.xlim([7, 16.5])
    plt.ylim([9, 16.5])
    plt.savefig('%s%s/%s_cell-cycle_%s_update1.pdf' % (output_dir, batch, batch, sample))
    plt.close()

    df.loc[len(df.index)] = [sample, identity, len(data_G1)/len(data), len(data_G1S)/len(data), len(data_S)/len(data), len(data_G2M)/len(data), len(data_G2MG1)/len(data)]

    plt.subplots(figsize=(9, 7))
    plt.hist([data_G1['hoechst'], data_S['hoechst'], data_G2M['hoechst']], weights=[np.ones(len(data_G1)) / len(data_G1), np.ones(len(data_S)) / len(data_S), np.ones(len(data_G2M)) / len(data_G2M)],
             range=[7000, 42000], bins=30, color=['#bc4d4a','#669daa', '#dd933c'], edgecolor='w')
    plt.savefig('%s/%s/%s_cell-cycle_hoechst_%s_update1.pdf' % (output_dir, batch, batch, sample))
    plt.close()

df.to_csv('%s/%s/%s_cell-cycle_update1.txt' % (data_dir, batch, batch), index=False, sep='\t')

