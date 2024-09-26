# Fig 1f

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import utilities as uti
import seaborn as sns

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
data_dir1 = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/processed/"
# output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/figures/"
batch = '5uM_24hr'
plate = 1
sample = 'XY165'
mode = 'all'  # only accepts 'neg', 'pos', 'all'

data_batch = pd.read_csv('%s/01_summary/%s.txt' % (data_dir, batch), na_values=['.'], sep='\t')
temp = data_batch[(data_batch['plate'] == plate) & (data_batch['sample'] == sample)]
hc = [temp['hoechst_cutoff_low'].tolist()[0], temp['hoechst_cutoff_high'].tolist()[0]]
cutoff = temp['log10_emiRFP670_cutoff'].tolist()[0]
red_cutoff = temp['red_cutoff'].tolist()[0]
green_cutoff = temp['green_cutoff'].tolist()[0]
data = pd.read_csv('%s/%s/%s_%s/txt/%s.txt' % (data_dir1, batch, batch, plate, sample), na_values=['.'], sep='\t')
data['log10_emiRFP670'] = np.log10(data['emiRFP670'])
data['log2_H2-3'] = np.log2(data['H2-3'])
data['log2_AzaleaB5'] = np.log2(data['AzaleaB5'])
data_filtered = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])].copy().reset_index(drop=True)
data_neg = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] < cutoff)].copy().reset_index(drop=True)
data_pos = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] >= cutoff)].copy().reset_index(drop=True)
if mode == 'pos':
    data_G1, data_G1S, data_S, data_G2M, data_G2MG1 = uti.get_cc(data_pos, red_cutoff, green_cutoff)
elif mode == 'neg':
    data_G1, data_G1S, data_S, data_G2M, data_G2MG1 = uti.get_cc(data_neg, red_cutoff, green_cutoff)
else:
    data_G1, data_G1S, data_S, data_G2M, data_G2MG1 = uti.get_cc(data_filtered, red_cutoff, green_cutoff)

fig, ax = plt.subplots(figsize=(9, 7))
# fig.subplots_adjust(right=0.8)
sns.scatterplot(data=data, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#cccccc')
sns.scatterplot(data=data_G1, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#bc4d4a')
sns.scatterplot(data=data_G1S, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#3d5a80')
sns.scatterplot(data=data_S, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#669daa')
sns.scatterplot(data=data_G2M, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#dd933c')
sns.scatterplot(data=data_G2MG1, x='log2_AzaleaB5', y='log2_H2-3', s=8, alpha=1, color='#d2b48c')
plt.xlim([7, 16.5])
plt.ylim([9, 16.5])
plt.show()

plt.subplots(figsize=(9, 7))
plt.hist([data_G1['hoechst'], data_S['hoechst'], data_G2M['hoechst']], weights=[np.ones(len(data_G1)) / len(data_G1), np.ones(len(data_S)) / len(data_S), np.ones(len(data_G2M)) / len(data_G2M)],
         range=[7000, 42000], bins=30, color=['#bc4d4a', '#669daa', '#dd933c'], edgecolor='w')
plt.show()

print("DONE!")