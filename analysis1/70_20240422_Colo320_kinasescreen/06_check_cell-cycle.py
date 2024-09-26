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

ref = pd.read_excel('%s/kinase_screen.xlsx' % master_folder, na_values=['.'])
print(ref.head())
hc = [5000, 40000]
cutoff = 2.95
n_fov = 9

batch = '5uM_24hr'
plate = 1
samples = ['XY02', 'XY04', 'XY06', 'XY08', 'XY10', 'XY68', 'XY12', 'XY14', 'XY16', 'XY18', 'XY20']

for i in range(len(samples)):
    print(i)
    sample = samples[i]
    data = pd.read_csv('%s/%s/%s_%s/txt/%s.txt' % (output_dir, batch, batch, plate, sample), na_values=['.'], sep='\t')
    data['log10_emiRFP670'] = np.log10(data['emiRFP670'])
    data['log2_H2-3'] = np.log2(data['H2-3'])
    data['log2_AzaleaB5'] = np.log2(data['AzaleaB5'])

    green_cutoff = 10.8
    red_cutoff = 9
    data_flt = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])]
    G2M = data_flt[(data_flt['log2_H2-3'] >= green_cutoff) & (data_flt['log2_AzaleaB5'] >= red_cutoff)]
    G1 = data_flt[(data_flt['log2_H2-3'] < green_cutoff) & (data_flt['log2_AzaleaB5'] >= red_cutoff)]
    G1S = data_flt[(data_flt['log2_H2-3'] < green_cutoff) & (data_flt['log2_AzaleaB5'] < red_cutoff)]
    S = data_flt[(data_flt['log2_H2-3'] >= green_cutoff) & (data_flt['log2_AzaleaB5'] < red_cutoff)]
    print(len(G1)/len(data_flt))
    print(len(G1S)/len(data_flt))
    print(len(S)/len(data_flt))
    print(len(G2M)/len(data_flt))

    fig, ax = plt.subplots(figsize=(8, 8))
    fig.subplots_adjust(right=0.8)
    sns.scatterplot(data=data_flt, x='log2_AzaleaB5', y='log2_H2-3', s=4, alpha=0.5)
    plt.axvline(x=red_cutoff, linestyle='--', color='r')
    plt.axhline(y=green_cutoff, linestyle='--', color='r')
    plt.xlim([7, 16.5])
    plt.ylim([9, 16.5])
    """if not os.path.exists("%s%s/%s_%s/hoechst_hist/" % (output_dir, batch, batch, plate)):
        os.makedirs("%s%s/%s_%s/hoechst_hist/" % (output_dir, batch, batch, plate))
    plt.savefig('%s%s/%s_%s/hoechst_hist/%s.pdf' % (output_dir, batch, batch, plate, sample))"""
    plt.show()



