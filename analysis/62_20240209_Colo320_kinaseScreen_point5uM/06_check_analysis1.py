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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240209_analysis_Colo320_kinaseScreen_point5uM/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

folders = ['R2P2', 'R2P3']
rep = 2
cell = 'HF+DFH'

target = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % data_dir, na_values=['.'])

data = pd.DataFrame()
for folder in folders:
    df = pd.read_csv('%s/%s/%s.txt' % (data_dir, folder, folder), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0)
data_screen = data[data['cell'] == cell].copy().reset_index(drop=True)
data_screen['HSR/DM'] = (data_screen['n_neg']+0.01)/(data_screen['n_pos']+0.01)
ctrl = np.mean(data_screen[data_screen['treatment'] == 'DMSO']['HSR/DM'])
data_screen['log2FC_HSR-DM'] = np.log2(data_screen['HSR/DM']/ctrl)

plt.subplots(figsize=(9, 6))
feature = 'log2FC_HSR-DM'
data_flt = data_screen[(data_screen['log2FC_HSR-DM'] < -0.7) | (data_screen['log2FC_HSR-DM'] > 0.5) |
                       (data_screen['n_total'] < 350)].copy().reset_index(drop=True)
print(data_flt['treatment'])
sns.scatterplot(data=data_screen, x=feature, y='n_total', alpha=0.5, s=7)
sns.scatterplot(data=data_screen[data_screen['treatment'] == 'DMSO'], x=feature, y='n_total', color='r', alpha=0.5, s=7)
for i in range(len(data_flt)):
    plt.text(x=data_flt['log2FC_HSR-DM'][i]+0.005, y=data_flt['n_total'][i]+0.005, s=data_flt['treatment'][i], size=5, color=(0/255, 191/255, 255/255))
plt.savefig('%s/R%s_%s_%s.pdf' % (output_dir, rep, cell, feature))
plt.show()

df = pd.DataFrame(columns=['screen', 'rep', 'cell', 'treatment', 'target', 'category'])
df['screen'] = [data['screen'].tolist()[0]] * len(data_flt)
df['rep'] = [rep] * len(data_flt)
df['cell'] = [cell] * len(data_flt)
df['treatment'] = data_flt['treatment']
df['target'] = [target[target['compound name'] == data_flt['treatment'][i]]['MedChemExpressTargets'].tolist()[0] for i in range(len(data_flt))]
category_lst = []
for i in range(len(data_flt)):
    if (data_flt['log2FC_HSR-DM'][i] < 0) & (data_flt['n_total'][i] < 250):
        category_lst.append('HSR_survival')
    elif data_flt['log2FC_HSR-DM'][i] < 0:
        category_lst.append('HSR')
    elif (data_flt['log2FC_HSR-DM'][i] > 0) & (data_flt['n_total'][i] < 250):
        category_lst.append('DM_survival')
    elif data_flt['log2FC_HSR-DM'][i] > 0:
        category_lst.append('DM')
    else:
        category_lst.append('NA')
df['category'] = category_lst
df.to_csv('%s/R%s_%s_hit.txt' % (output_dir, rep, cell), index=False, sep='\t')




