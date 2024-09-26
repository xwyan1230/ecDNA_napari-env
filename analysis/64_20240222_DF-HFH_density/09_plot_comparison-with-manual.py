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

folders = ['48hr_384well']
test_lst = ['XY14', 'XY15', 'XY14_manual', 'XY15_manual']
density_lst = [10000, 8000, 6000, 4000, 3000, 2000, 1000]
size_lst = [40, 35,30,25,20,15,10]
ctrl_colors = [(0.6, 0.81, 0.98), (0.5, 0.81, 0.98), (0.4, 0.81, 0.98), (0.3, 0.81, 0.98), (0.2, 0.81, 0.98), (0.1, 0.81, 0.98), (0, 0.81, 0.98)][::-1]
MK1775_colors = [(0.88, 0.6, 0.6), (0.88, 0.5, 0.5), (0.88, 0.4, 0.4), (0.88, 0.3, 0.3), (0.88, 0.2, 0.2), (0.88, 0.1, 0.1), (0.88, 0, 0)][::-1]
rainboo_colors = [(220/255, 20/255, 60/255), (1, 140/255, 0), (1, 215/255, 0), (154/255, 205/255, 50/255), (64/255, 224/255, 208/255),
                  (100/255, 149/255, 237/255), (123/255, 104/255, 238/255)]


data = pd.DataFrame()
for folder in folders:
    df = pd.read_csv('%s/%s/%s_compare_fov.txt' % (data_dir, folder, folder), na_values=['.'], sep='\t')
    df['treatment'] = [df['sample_name'][x].split('_')[0] for x in range(len(df))]
    df['HSR/DM'] = (df['n_pos']+0.01)/(df['n_neg']+0.01)
    data = pd.concat([data, df], axis=0)

plt.subplots(figsize=(9, 6))
feature = 'HSR/DM'
sns.scatterplot(data=data[data['sample'].isin(['XY14'])], x=feature, y='n_total', alpha=0.8, s=40, color=(0.88, 0.1, 0.1), label='MK1775_auto')
sns.scatterplot(data=data[data['sample'].isin(['XY15'])], x=feature, y='n_total', alpha=0.8, s=40, color=(0.88, 0.4, 0.4), label='DMSO_auto')
sns.scatterplot(data=data[data['sample'].isin(['XY14_manual'])], x=feature, y='n_total', alpha=0.8, s=40, color=(0.1, 0.81, 0.98), label='MK1775_manual')
sns.scatterplot(data=data[data['sample'].isin(['XY15_manual'])], x=feature, y='n_total', alpha=0.8, s=40, color=(0.4, 0.81, 0.98), label='DMSO_manual')
plt.xlim([0, 2.5])
plt.ylim([0, 1200])
plt.legend()
plt.savefig('%s/%s_comparison_auto-vs-manual.pdf' % (output_dir, folder))
plt.show()

plt.subplots(figsize=(9, 6))
feature = 'HSR/DM'
sns.scatterplot(data=data[data['sample'].isin(['XY14_update'])], x=feature, y='n_total', alpha=0.8, s=40, color=(0.88, 0.1, 0.1), label='MK1775_auto1')
sns.scatterplot(data=data[data['sample'].isin(['XY15_update'])], x=feature, y='n_total', alpha=0.8, s=40, color=(0.88, 0.4, 0.4), label='DMSO_auto1')
sns.scatterplot(data=data[data['sample'].isin(['XY14_manual'])], x=feature, y='n_total', alpha=0.8, s=40, color=(0.1, 0.81, 0.98), label='MK1775_manual')
sns.scatterplot(data=data[data['sample'].isin(['XY15_manual'])], x=feature, y='n_total', alpha=0.8, s=40, color=(0.4, 0.81, 0.98), label='DMSO_manual')
plt.xlim([0, 2.5])
plt.ylim([0, 1200])
plt.legend()
plt.savefig('%s/%s_comparison_auto1-vs-manual.pdf' % (output_dir, folder))
plt.show()

plt.subplots(figsize=(9, 6))
feature = 'HSR/DM'
sns.scatterplot(data=data[data['sample'].isin(['XY14_update1'])], x=feature, y='n_total', alpha=0.8, s=40, color=(0.88, 0.1, 0.1), label='MK1775_auto2')
sns.scatterplot(data=data[data['sample'].isin(['XY15_update1'])], x=feature, y='n_total', alpha=0.8, s=40, color=(0.88, 0.4, 0.4), label='DMSO_auto2')
sns.scatterplot(data=data[data['sample'].isin(['XY14_manual'])], x=feature, y='n_total', alpha=0.8, s=40, color=(0.1, 0.81, 0.98), label='MK1775_manual')
sns.scatterplot(data=data[data['sample'].isin(['XY15_manual'])], x=feature, y='n_total', alpha=0.8, s=40, color=(0.4, 0.81, 0.98), label='DMSO_manual')
plt.xlim([0, 2.5])
plt.ylim([0, 1200])
plt.legend()
plt.savefig('%s/%s_comparison_auto2-vs-manual.pdf' % (output_dir, folder))
plt.show()






