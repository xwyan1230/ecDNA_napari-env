import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import skimage.io as skio
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
data_dir1 = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

# ref = pd.read_excel('%s/kinase_screen.xlsx' % master_folder, na_values=['.'])
# print(ref.head())
hc = [5000, 40000]
cutoff = 2.95
n_fov = 9

batch = '48hr_density'

skip = pd.read_csv('%s/skip.txt' % data_dir1, na_values=['.'], sep='\t')

samples = ['XY0%s' % (x+1) for x in range(9)] + ['XY%s' % (x+10) for x in range(201)]
wells = ['D%s' % (x + 3) for x in range(20)] + ['E%s' % (x + 3) for x in range(20)][::-1] + \
            ['F%s' % (x + 3) for x in range(20)] + ['G%s' % (x + 3) for x in range(20)][::-1] + \
            ['H%s' % (x + 3) for x in range(20)] + ['I%s' % (x + 3) for x in range(20)][::-1] + \
            ['J%s' % (x + 3) for x in range(20)] + ['K%s' % (x + 3) for x in range(20)][::-1] + \
            ['L%s' % (x + 3) for x in range(20)] + ['M%s' % (x + 3) for x in range(20)][::-1] + \
            ['N%s' % (x + 3) for x in range(10)]
densities = ['1k'] * 20 + ['2k'] * 20 + ['3k'] * 20 + ['4k'] * 20 + ['6k'] * 20 + ['8k'] * 20 + \
            ['10k'] * 20 + ['12k'] * 20 + ['14k'] * 20 + ['16k'] * 20 + ['ctrl'] * 10


def img_crop(img, i):
    if i == 0:
        out = img[0:1440, 0:1920]
    elif (i == 1) | (i == 2):
        out = img[0:1440, 576:1920]
    elif (i == 3) | (i == 6):
        out = img[432:1440, 0:1920]
    elif (i == 4) | (i == 5):
        out = img[432:1440, 0:1344]
    elif (i == 7) | (i == 8):
        out = img[432:1440, 576:1920]
    return out


def img_factor(i):
    if i == 0:
        out = 1
    elif (i == 1) | (i == 2) | (i == 3) | (i == 6):
        out = 0.7
    elif (i == 4) | (i == 5) | (i == 7) | (i == 8):
        out = 0.49
    return out


df_fov = pd.DataFrame(columns=['screen', 'sample', 'well', 'density', 'hoechst_cutoff', 'log10_emiRFP670_cutoff',
                               'fov', 'hoechst_mean', 'n_total', 'n_filtered', 'n_neg', 'n_pos', 'per_neg', 'per_pos', 'fov_hoechst'])
df = pd.DataFrame(columns=['screen', 'sample', 'well', 'density', 'hoechst_cutoff', 'log10_emiRFP670_cutoff',
                   'hoechst_mean', 'n_total', 'n_filtered', 'n_neg', 'n_pos', 'per_neg', 'per_pos', 'fov_hoechst'])

for i in range(len(samples)):
    print(i)
    sample = samples[i]
    well = wells[i]
    density = densities[i]
    data = pd.read_csv('%s/%s/%s/txt/%s.txt' % (output_dir, batch, batch, sample), na_values=['.'], sep='\t')
    data['log10_emiRFP670'] = np.log10(data['emiRFP670'])

    fig, ax = plt.subplots(figsize=(9, 6))
    fig.subplots_adjust(right=0.8)
    sns.histplot(data=data, x='hoechst')
    plt.axvline(hc[0], 0, 1000, c='r')
    plt.axvline(hc[1], 0, 1000, c='r')
    if not os.path.exists("%s%s/%s/hoechst_hist/" % (output_dir, batch, batch)):
        os.makedirs("%s%s/%s/hoechst_hist/" % (output_dir, batch, batch))
    plt.savefig('%s%s/%s/hoechst_hist/%s.pdf' % (output_dir, batch, batch, sample))
    plt.close()

    fig, ax = plt.subplots(figsize=(9, 6))
    sns.histplot(data=data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])], x='log10_emiRFP670', bins=20)  # annot=True
    plt.xlim([2, 6])
    plt.axvline(cutoff, 0, 1000, c='r')
    if not os.path.exists("%s%s/%s/emiRFP670_hist/" % (output_dir, batch, batch)):
        os.makedirs("%s%s/%s/emiRFP670_hist/" % (output_dir, batch, batch))
    plt.savefig("%s%s/%s/emiRFP670_hist/%s.pdf" % (output_dir, batch, batch, sample))
    plt.close()

    for j in range(n_fov):
        if '%s_%s_%s' % (batch, sample, j + 1) not in skip['samples'].tolist():
            file_name = 'Image_%s_0000%s' % (sample, j + 1)
            img_hoechst = img_crop(skio.imread("%s%s/%s/%s/%s_CH1.tif" % (data_dir1, batch, batch, sample, file_name), plugin="tifffile")[:, :, 2], j)
            fov_hoechst = np.sum(img_hoechst) / img_factor(j)
            temp = data[data['fov'] == j].copy().reset_index(drop=True)
            hoechst_mean = np.mean(temp['hoechst'])
            n_total = len(temp)
            n_filtered = len(temp[(temp['hoechst'] > hc[0]) & (temp['hoechst'] < hc[1])])
            n_neg = len(temp[(temp['hoechst'] > hc[0]) & (temp['hoechst'] < hc[1]) & (temp['log10_emiRFP670'] < cutoff)])
            n_pos = len(temp[(temp['hoechst'] > hc[0]) & (temp['hoechst'] < hc[1]) & (temp['log10_emiRFP670'] >= cutoff)])
            per_neg = n_neg / (n_filtered+0.01)
            per_pos = n_pos / (n_filtered+0.01)
            df_fov.loc[len(df_fov.index)] = [data['screen'][0], sample, well, density, hc, cutoff,
                                             j + 1, hoechst_mean, n_total, n_filtered, n_neg, n_pos, per_neg, per_pos, fov_hoechst]

    hoechst_mean = np.mean(data['hoechst'])
    n_total = len(data)
    n_filtered = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])])
    n_neg = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] < cutoff)])
    n_pos = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] >= cutoff)])
    per_neg = n_neg/(n_filtered+0.01)
    per_pos = n_pos/(n_filtered+0.01)
    fov_hoechst_mean = np.mean(df_fov[df_fov['sample'] == sample]['fov_hoechst'])

    df.loc[len(df.index)] = [data['screen'][0], sample, well, density, hc, cutoff, hoechst_mean,
                             n_total, n_filtered, n_neg, n_pos, per_neg, per_pos, fov_hoechst_mean]

df.to_csv('%s/%s/%s.txt' % (output_dir, batch, batch), index=False, sep='\t')
df_fov.to_csv('%s/%s/%s_fov.txt' % (output_dir, batch, batch), index=False, sep='\t')



