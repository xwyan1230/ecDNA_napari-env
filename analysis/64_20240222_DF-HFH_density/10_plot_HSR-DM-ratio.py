import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
from shared.sinaplot import sinaplot
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import scipy.stats as stats
import shared.dataframe as dat
import shared.math as mat
import math
import seaborn as sns
from scipy.stats import iqr
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240222_analysis_DF-HFH_density/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

folders = ['72hr_384well']
test_lst = [['XY0%s' % (x+3) for x in range(7)] + ['XY%s' % (x+10) for x in range(3)],
             ['XY%s' % (x+13) for x in range(10)], ['XY%s' % (x+23) for x in range(10)],
             ['XY%s' % (x+33) for x in range(10)], ['XY%s' % (x+43) for x in range(10)],
             ['XY%s' % (x+53) for x in range(10)], ['XY%s' % (x+63) for x in range(10)]]
density_lst = [10000, 8000, 6000, 4000, 3000, 2000, 1000]
size_lst = [40, 35,30,25,20,15,10]
ctrl_colors = [(0.6, 0.81, 0.98), (0.5, 0.81, 0.98), (0.4, 0.81, 0.98), (0.3, 0.81, 0.98), (0.2, 0.81, 0.98), (0.1, 0.81, 0.98), (0, 0.81, 0.98)][::-1]
MK1775_colors = [(0.88, 0.6, 0.6), (0.88, 0.5, 0.5), (0.88, 0.4, 0.4), (0.88, 0.3, 0.3), (0.88, 0.2, 0.2), (0.88, 0.1, 0.1), (0.88, 0, 0)][::-1]
rainboo_colors = [(220/255, 20/255, 60/255), (1, 140/255, 0), (1, 215/255, 0), (154/255, 205/255, 50/255), (64/255, 224/255, 208/255),
                  (100/255, 149/255, 237/255), (123/255, 104/255, 238/255)]


data = pd.DataFrame()
for folder in folders:
    df = pd.read_csv('%s/%s/%s_update1_hc.txt' % (data_dir, folder, folder), na_values=['.'], sep='\t')
    df['treatment'] = [df['sample_name'][x].split('_')[0] for x in range(len(df))]
    df['HSR/DM'] = (df['n_pos']+0.01)/(df['n_neg']+0.01)
    data = pd.concat([data, df], axis=0)

plt.subplots(figsize=(9, 6))
feature = 'HSR/DM'
for i in range(len(test_lst)):
    test = test_lst[i]
    data_screen = data[data['sample'].isin(test)].copy().reset_index(drop=True)
    ctrl = np.mean(data_screen[data_screen['treatment'] == 'DMSO']['HSR/DM'])
    data_screen['log2FC_HSR-DM'] = np.log2(data_screen['HSR/DM']/ctrl)

    sns.scatterplot(data=data_screen, x=feature, y='n_total', alpha=0.8, s=size_lst[i], color=MK1775_colors[i])
    sns.scatterplot(data=data_screen[data_screen['treatment'] == 'DMSO'], x=feature, y='n_total', color=ctrl_colors[i], alpha=0.8, s=size_lst[i])
plt.xlim([0, 3])
plt.ylim([0, 5000])
# plt.savefig('%s/%s_MK1775_update.pdf' % (output_dir, folder))
plt.show()

plt.subplots(figsize=(9, 6))
for i in range(len(test_lst)):
    test = test_lst[i]
    data_screen = data[data['sample'].isin(test)].copy().reset_index(drop=True)
    ctrl = np.mean(data_screen[data_screen['treatment'] == 'DMSO']['HSR/DM'])
    data_screen['log2FC_HSR-DM'] = np.log2(data_screen['HSR/DM']/ctrl)

    sns.scatterplot(data=data_screen, x=feature, y='n_total', alpha=0.8, s=size_lst[i], color=rainboo_colors[i])
plt.xlim([0, 3])
plt.ylim([0, 5000])
# plt.savefig('%s/%s_density_update.pdf' % (output_dir, folder))
plt.show()






