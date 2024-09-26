import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240222_analysis_DF-HFH_density/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sprocessed/" % master_folder

hc = [5000, 40000]
cutoff = 3.1
n_fov = 9

folder = '48hr_384well'
sample_names = ['MK1775_4x8000-rep1', 'MK1775_4x8000-rep1_update', 'MK1775_4x8000-rep1_update1', 'MK1775_4x8000-rep1_manual',
                'DMSO_4x8000-rep8', 'DMSO_4x8000-rep8_update',  'DMSO_4x8000-rep8_update1', 'DMSO_4x8000-rep8_manual']

samples = ['XY14', 'XY14_update', 'XY14_update1', 'XY14_manual', 'XY15', 'XY15_update', 'XY15_update1', 'XY15_manual']

df_fov = pd.DataFrame(columns=['sample', 'sample_name', 'hoechst_cutoff', 'log10_emiRFP670_cutoff',
                               'fov', 'hoechst_mean', 'n_total', 'n_filtered', 'n_neg', 'n_pos', 'per_neg', 'per_pos'])
df = pd.DataFrame(columns=['sample', 'sample_name', 'hoechst_cutoff', 'log10_emiRFP670_cutoff',
                   'hoechst_mean', 'n_total', 'n_filtered', 'n_neg', 'n_pos', 'per_neg', 'per_pos'])

for i in range(len(samples)):
    print(i)
    sample = samples[i]
    sample_name = sample_names[i]
    data = pd.read_csv('%s/%s/txt/%s.txt' % (data_dir, folder, sample_name), na_values=['.'], sep='\t')
    data['log10_emiRFP670'] = np.log10(data['emiRFP670'])

    """fig, ax = plt.subplots(figsize=(9, 6))
    fig.subplots_adjust(right=0.8)
    sns.histplot(data=data, x='hoechst')
    plt.axvline(hc, 0, 1000, c='r')
    if not os.path.exists("%s%s/hoechst_hist/" % (output_dir, folder)):
        os.makedirs("%s%s/hoechst_hist/" % (output_dir, folder))
    plt.savefig('%s/%s/hoechst_hist/%s.pdf' % (output_dir, folder, sample))
    plt.close()

    fig, ax = plt.subplots(figsize=(9, 6))
    sns.histplot(data=data[data['hoechst'] > hc], x='log10_emiRFP670', bins=20)  # annot=True
    plt.xlim([2, 6])
    plt.axvline(cutoff, 0, 1000, c='r')
    if not os.path.exists("%s%s/emiRFP670_hist/" % (output_dir, folder)):
        os.makedirs("%s%s/emiRFP670_hist/" % (output_dir, folder))
    plt.savefig('%s/%s/emiRFP670_hist/%s.pdf' % (output_dir, folder, sample))
    plt.close()"""

    hoechst_mean = np.mean(data['hoechst'])
    n_total = len(data)
    n_filtered = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])])
    n_neg = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] < cutoff)])
    n_pos = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] >= cutoff)])
    per_neg = n_neg/(n_filtered+0.01)
    per_pos = n_pos/(n_filtered+0.01)

    df.loc[len(df.index)] = [sample,
                             sample_name, hc, cutoff, hoechst_mean, n_total, n_filtered, n_neg, n_pos, per_neg, per_pos]

    for j in range(n_fov):
        temp = data[data['fov'] == j].copy().reset_index(drop=True)
        hoechst_mean = np.mean(temp['hoechst'])
        n_total = len(temp)
        n_filtered = len(temp[(temp['hoechst'] > hc[0]) & (temp['hoechst'] < hc[1])])
        n_neg = len(temp[(temp['hoechst'] > hc[0]) & (temp['hoechst'] < hc[1]) & (temp['log10_emiRFP670'] < cutoff)])
        n_pos = len(temp[(temp['hoechst'] > hc[0]) & (temp['hoechst'] < hc[1]) & (temp['log10_emiRFP670'] >= cutoff)])
        per_neg = n_neg / (n_filtered+0.01)
        per_pos = n_pos / (n_filtered+0.01)
        df_fov.loc[len(df_fov.index)] = [sample, sample_name, hc, cutoff, j+1, hoechst_mean, n_total, n_filtered,
                                         n_neg, n_pos, per_neg, per_pos]

df.to_csv('%s/%s/%s_compare.txt' % (output_dir, folder, folder), index=False, sep='\t')
df_fov.to_csv('%s/%s/%s_compare_fov.txt' % (output_dir, folder, folder), index=False, sep='\t')



