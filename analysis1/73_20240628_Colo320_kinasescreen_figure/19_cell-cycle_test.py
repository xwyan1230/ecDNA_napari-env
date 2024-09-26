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

folder = '5uM_48hr'
batch = '5uM_48hr_1'
if int(batch[-1]) == 1:
    samples = ['XY02', 'XY04', 'XY06', 'XY08', 'XY10', 'XY12', 'XY14', 'XY16', 'XY18', 'XY20']
    samples = ['XY156']
    samples = ['XY0%s' % (x+1) for x in range(9)] + ['XY%s' % (x+10) for x in range(191)]
elif int(batch[-1]) == 2:
    samples = ['XY%s' % (2*x+111) for x in range(5)] + ['XY%s' % (2*x+101) for x in range(5)]
    samples = ['XY178', 'XY185']
    # samples = ['XY70', 'XY80', 'XY90', 'XY150', 'XY170']
    # samples = ['XY02', 'XY04', 'XY06', 'XY08', 'XY10', 'XY12', 'XY14', 'XY16', 'XY18', 'XY20']
else:
    samples = ['XY01', 'XY02']+['XY%s' % (2*x+201) for x in range(5)] + ['XY%s' % (2*x+202) for x in range(5)]
    samples = ['XY45']
    # samples = ['XY0%s' % (x+1) for x in range(9)] + ['XY%s' % (x+10) for x in range(5)]

identities = ['neg'] * 5 + ['pos'] * 5 + ['mix'] * 40
# colors = ['#808080', '#bc4d4a'] * 5 + []


def get_ellipse(x, y, a, b, A, data, data_x, data_y):
    data_flt = data[(((data[data_x]-x)*math.cos(A) + (data[data_y]-y)*math.sin(A))**2)/(a**2)+(((data[data_x]-x)*math.sin(A) + (data[data_y]-y)*math.cos(A))**2)/(b**2) <=1].copy().reset_index(drop=True)
    return data_flt


def get_ellipse1(x, y, a, b, A, data, data_x, data_y):
    c = math.sqrt(a**2-b**2)
    data_flt = data[np.sqrt((data[data_x]-x+c*math.cos(A))**2+(data[data_y]-y+c*math.sin(A))**2)+np.sqrt((data[data_x]-x-c*math.cos(A))**2+(data[data_y]-y-c*math.sin(A))**2)-2*a<=0].copy().reset_index(drop=True)
    return data_flt


df = pd.DataFrame(columns=['sample', 'identity', 'G1', 'S', 'G2M'])
for i in range(len(samples)):
    sample = samples[i]
    print(sample)
    # identity = identities[i]
    data = pd.read_csv('%s/%s/%s/txt/%s.txt' % (data_dir, folder, batch, sample), na_values=['.'], sep='\t')
    data['log10_emiRFP670'] = np.log10(data['emiRFP670'])
    data['log2_H2-3'] = np.log2(data['H2-3'])
    data['log2_AzaleaB5'] = np.log2(data['AzaleaB5'])
    axis = np.median(data[data['log2_AzaleaB5']<9]['log2_AzaleaB5'])-0.1
    axis=8.3
    print(axis)

    data = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])].reset_index(drop=True)
    data_G1 = get_ellipse1(12.3, 10.1, 3.3, 0.7, 0, data, 'log2_AzaleaB5', 'log2_H2-3')  # 24hr_density
    # data_G1 = get_ellipse1(13, 10.1, 4, 0.7, math.radians(6), data, 'log2_AzaleaB5', 'log2_H2-3')  # 24hr_density
    # data_G1 = get_ellipse1(12.6, 10.1, 3.2, 0.7, 0, data, 'log2_AzaleaB5', 'log2_H2-3')
    data_G1['group'] = ['G1']* len(data_G1)
    data_S = get_ellipse1(axis, 13.2, 2.5, 0.7, math.radians(85), data, 'log2_AzaleaB5', 'log2_H2-3')  # 24hr_density
    # data_S = get_ellipse1(8.9, 13.2, 2.5, 0.7, math.radians(85), data, 'log2_AzaleaB5', 'log2_H2-3')
    data_S['group'] = ['S'] * len(data_S)
    data_G2M = get_ellipse1(12+axis-8.3, 14, 3.5, 1.3, math.radians(30), data, 'log2_AzaleaB5', 'log2_H2-3')  # 24hr_density
    # data_G2M = get_ellipse1(12.6, 14, 3.5, 1.3, math.radians(30), data, 'log2_AzaleaB5', 'log2_H2-3')
    data_G2M['group'] = ['G2M'] * len(data_G2M)
    data_all = pd.concat([data_G1, data_S, data_G2M], axis=0).reset_index(drop=True)
    # df.loc[len(df.index)] = [sample, identity, len(data_G1)/len(data_all), len(data_S)/len(data_all), len(data_G2M)/len(data_all)]

    fig, ax = plt.subplots(figsize=(9, 7))
    # fig.subplots_adjust(right=0.8)
    sns.scatterplot(data=data, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#cccccc') ##4d4e4e
    sns.scatterplot(data=data_G1, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#bc4d4a')
    sns.scatterplot(data=data_S, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#669daa')
    sns.scatterplot(data=data_G2M, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#dd933c')
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

