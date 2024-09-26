import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import skimage.io as skio
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240222_analysis_DF-HFH_density/"
data_dir = "%sprocessed/" % master_folder
data_dir1 = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

hc = [5000, 40000]
# hc = [0, 66000]
cutoff = 3.1
n_fov = 9


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


folder = '24hr_96well'
sample_names = ['DF-Ctrl', 'HFH-Ctrl'] + \
               ['DMSO_4x10000-rep%s' % (x+1) for x in range(8)] + \
               ['MK1775_4x10000-rep%s' % (x+1) for x in range(2)] + \
               ['MK1775_4x8000-rep%s' % (x+1) for x in range(2)][::-1] + \
               ['DMSO_4x8000-rep%s' % (x+1) for x in range(8)][::-1] + \
               ['DMSO_4x6000-rep%s' % (x+1) for x in range(8)] + \
               ['MK1775_4x6000-rep%s' % (x+1) for x in range(2)] + \
               ['MK1775_4x4000-rep%s' % (x+1) for x in range(2)][::-1] + \
               ['DMSO_4x4000-rep%s' % (x+1) for x in range(8)][::-1] + \
               ['DMSO_4x3000-rep%s' % (x+1) for x in range(8)] + \
               ['MK1775_4x3000-rep%s' % (x+1) for x in range(2)] + \
               ['MK1775_4x2000-rep%s' % (x+1) for x in range(2)][::-1] + \
               ['DMSO_4x2000-rep%s' % (x+1) for x in range(8)][::-1] + \
               ['DMSO_4x1000-rep%s' % (x+1) for x in range(8)] + \
               ['MK1775_4x1000-rep%s' % (x+1) for x in range(2)]

samples = ['XY0%s' % (x+1) for x in range(9)] + ['XY%s' % (x+10) for x in range(63)]

df_fov = pd.DataFrame(columns=['exp', 'group', 'timepoint', 'format', 'sample', 'sample_name', 'hoechst_cutoff', 'log10_emiRFP670_cutoff',
                               'fov', 'hoechst_mean', 'n_total', 'n_filtered', 'n_neg', 'n_pos', 'per_neg', 'per_pos', 'fov_hoechst'])
df = pd.DataFrame(columns=['exp', 'group', 'timepoint', 'format', 'sample', 'sample_name', 'hoechst_cutoff', 'log10_emiRFP670_cutoff',
                   'hoechst_mean', 'n_total', 'n_filtered', 'n_neg', 'n_pos', 'per_neg', 'per_pos', 'fov_hoechst'])

for i in range(len(samples)):
    print(i)
    sample = samples[i]
    sample_name = sample_names[i]
    data = pd.read_csv('%s/%s/txt/%s_update1.txt' % (data_dir, folder, sample_name), na_values=['.'], sep='\t')
    data['log10_emiRFP670'] = np.log10(data['emiRFP670'])

    fig, ax = plt.subplots(figsize=(9, 6))
    fig.subplots_adjust(right=0.8)
    sns.histplot(data=data, x='hoechst')
    plt.axvline(hc[0], 0, 1000, c='r')
    plt.axvline(hc[1], 0, 1000, c='r')
    if not os.path.exists("%s%s/hoechst_hist_update1_hc/" % (output_dir, folder)):
        os.makedirs("%s%s/hoechst_hist_update1_hc/" % (output_dir, folder))
    plt.savefig('%s/%s/hoechst_hist_update1_hc/%s.pdf' % (output_dir, folder, sample))
    plt.close()

    fig, ax = plt.subplots(figsize=(9, 6))
    sns.histplot(data=data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])], x='log10_emiRFP670', bins=20)  # annot=True
    plt.xlim([2, 6])
    plt.axvline(cutoff, 0, 1000, c='r')
    if not os.path.exists("%s%s/emiRFP670_hist_update1_hc/" % (output_dir, folder)):
        os.makedirs("%s%s/emiRFP670_hist_update1_hc/" % (output_dir, folder))
    plt.savefig('%s/%s/emiRFP670_hist_update1_hc/%s.pdf' % (output_dir, folder, sample))
    plt.close()

    for j in range(n_fov):
        file_name = 'Image_%s_0000%s' % (sample, j + 1)
        img_hoechst = img_crop(skio.imread("%s%s/%s/%s_CH1.tif" % (data_dir1, folder, sample, file_name), plugin="tifffile")[:, :, 2], j)
        fov_hoechst = np.sum(img_hoechst) / img_factor(j)
        temp = data[data['fov'] == j].copy().reset_index(drop=True)
        hoechst_mean = np.mean(temp['hoechst'])
        n_total = len(temp)
        n_filtered = len(temp[(temp['hoechst'] > hc[0]) & (temp['hoechst'] < hc[1])])
        n_neg = len(temp[(temp['hoechst'] > hc[0]) & (temp['hoechst'] < hc[1]) & (temp['log10_emiRFP670'] < cutoff)])
        n_pos = len(temp[(temp['hoechst'] > hc[0]) & (temp['hoechst'] < hc[1]) & (temp['log10_emiRFP670'] >= cutoff)])
        per_neg = n_neg / (n_filtered+0.01)
        per_pos = n_pos / (n_filtered+0.01)
        df_fov.loc[len(df_fov.index)] = [data['exp'][0], data['group'][0], data['timepoint'][0], data['format'][0],
                                         sample, sample_name, hc, cutoff, j+1, hoechst_mean, n_total, n_filtered,
                                         n_neg, n_pos, per_neg, per_pos, fov_hoechst]

    hoechst_mean = np.mean(data['hoechst'])
    n_total = len(data)
    n_filtered = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])])
    n_neg = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] < cutoff)])
    n_pos = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] >= cutoff)])
    per_neg = n_neg/(n_filtered+0.01)
    per_pos = n_pos/(n_filtered+0.01)
    fov_hoechst_mean = np.mean(df_fov[df_fov['sample'] == sample]['fov_hoechst'])

    df.loc[len(df.index)] = [data['exp'][0], data['group'][0], data['timepoint'][0], data['format'][0], sample,
                             sample_name, hc, cutoff, hoechst_mean, n_total, n_filtered, n_neg, n_pos, per_neg, per_pos, fov_hoechst_mean]

df.to_csv('%s/%s/%s_update1_hc.txt' % (output_dir, folder, folder), index=False, sep='\t')
df_fov.to_csv('%s/%s/%s_update1_hc_fov.txt' % (output_dir, folder, folder), index=False, sep='\t')



