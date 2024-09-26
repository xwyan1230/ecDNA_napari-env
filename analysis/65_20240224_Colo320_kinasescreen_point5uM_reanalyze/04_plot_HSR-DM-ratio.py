import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from shared.sinaplot import sinaplot
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import random
import scipy.stats as stats
import shared.dataframe as dat
import shared.math as mat
import math
import seaborn as sns
from scipy.stats import iqr
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240209_analysis_Colo320_kinaseScreen_point5uM/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

folders = ['R2P1', 'R2P2']
rep = 2
cell = 'DF+HFH'
hc = [5000, 40000]
cutoff = 2.95

target = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % data_dir, na_values=['.'])
print(target.head())
ctrl_colors = [(0.6, 0.81, 0.98), (0.5, 0.81, 0.98), (0.4, 0.81, 0.98), (0.3, 0.81, 0.98), (0.2, 0.81, 0.98), (0.1, 0.81, 0.98), (0, 0.81, 0.98)][::-1]
MK1775_colors = [(0.88, 0.6, 0.6), (0.88, 0.5, 0.5), (0.88, 0.4, 0.4), (0.88, 0.3, 0.3), (0.88, 0.2, 0.2), (0.88, 0.1, 0.1), (0.88, 0, 0)][::-1]

data = pd.DataFrame()
for folder in folders:
    df = pd.read_csv('%s/%s/%s_update1_hc.txt' % (data_dir, folder, folder), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0)
data_screen = data[data['cell'] == cell].copy().reset_index(drop=True)
data_screen['HSR/DM'] = (data_screen['n_pos']+0.01)/(data_screen['n_neg']+0.01)

data_ctrl = data_screen[data_screen['treatment'] == 'DMSO'].copy().reset_index(drop=True)
data_ctrl_cells = pd.DataFrame()
for i in range(len(data_ctrl)):
    folder = 'R%sP%s' % (data_ctrl['rep'][i], data_ctrl['plate'][i])
    sample = data_ctrl['sample'][i]
    df = pd.read_csv('%s/%s/txt/%s_update1.txt' % (data_dir, folder, sample), na_values=['.'], sep='\t')
    df = df[(df['hoechst'] > hc[0]) & (df['hoechst'] < hc[1])].copy().reset_index(drop=True)
    data_ctrl_cells = pd.concat([data_ctrl_cells, df], axis=0)
data_ctrl_cells['log10_emiRFP670'] = np.log10(data_ctrl_cells['emiRFP670'])

data_ctrl_random = pd.DataFrame(columns=['screen', 'rep', 'cell', 'treatment', 'hoechst_cutoff', 'log10_emiRFP670_cutoff', 'n_filtered', 'n_neg', 'n_pos', 'per_neg', 'per_pos'])
for i in range(1000):
    print(i)
    n_total_lst = np.arange(1, 3600, 1)
    n_total_temp = random.choice(n_total_lst)
    df_temp = data_ctrl_cells.sample(n=n_total_temp)
    n_neg = len(df_temp[df_temp['log10_emiRFP670'] < cutoff])
    n_pos = len(df_temp[df_temp['log10_emiRFP670'] >= cutoff])
    per_neg = n_neg / (n_total_temp + 0.01)
    per_pos = n_pos / (n_total_temp + 0.01)
    data_ctrl_random.loc[len(data_ctrl_random.index)] = [data_ctrl_cells['screen'][0], rep, cell, 'DMSO', hc, cutoff,
                                                         n_total_temp, n_neg, n_pos, per_neg, per_pos]
data_ctrl_random['HSR/DM'] = (data_ctrl_random['n_pos']+0.01)/(data_ctrl_random['n_neg']+0.01)
data_ctrl_random.to_csv('%s/R%s_ctrl_random_update1_hc.txt' % (output_dir, rep), index=False, sep='\t')

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_ctrl_random, x='HSR/DM', y='n_filtered', alpha=0.5, s=30, color=ctrl_colors[6])
sns.scatterplot(data=data_screen, x='HSR/DM', y='n_filtered', alpha=0.5, s=30, color=MK1775_colors[0])
sns.scatterplot(data=data_screen[data_screen['treatment'] == 'DMSO'], x='HSR/DM', y='n_filtered', color=ctrl_colors[0], alpha=0.5, s=30)
std_ratio = np.std(data_screen[data_screen['treatment'] == 'DMSO']['HSR/DM'])
mean_ratio = np.mean(data_screen[data_screen['treatment'] == 'DMSO']['HSR/DM'])
plt.xlim([0, 3])
plt.ylim([0, 3600])
plt.legend()
# plt.savefig('%s/R%s_%s_n-neg_vs_n-pos.pdf' % (output_dir, rep, cell))
plt.show()








