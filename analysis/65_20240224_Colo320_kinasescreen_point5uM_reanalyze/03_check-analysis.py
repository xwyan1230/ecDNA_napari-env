import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
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

folders = ['R2P1', 'R2P2']
rep = 2
cell = 'DF+HFH'

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
# ctrl = np.mean(data_screen[data_screen['treatment'] == 'DMSO']['HSR/DM'])
# data_screen['log2FC_HSR-DM'] = np.log2(data_screen['HSR/DM']/ctrl)


def func(x, a):
    return a * x


plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_screen, x='n_neg', y='n_pos', alpha=0.5, s=30, color=MK1775_colors[0])
sns.scatterplot(data=data_screen[data_screen['treatment'] == 'DMSO'], x='n_neg', y='n_pos', color=ctrl_colors[0], alpha=0.5, s=30)
std_ratio = np.std(data_screen[data_screen['treatment'] == 'DMSO']['HSR/DM'])
mean_ratio = np.mean(data_screen[data_screen['treatment'] == 'DMSO']['HSR/DM'])
print(mean_ratio-3*std_ratio)
print(mean_ratio+3*std_ratio)
std_n_pos = np.std(data_screen[data_screen['treatment'] == 'DMSO']['n_pos'])
mean_n_pos = np.mean(data_screen[data_screen['treatment'] == 'DMSO']['n_pos'])
print(mean_n_pos-3*std_n_pos)
print(mean_n_pos+3*std_n_pos)
std_n_neg = np.std(data_screen[data_screen['treatment'] == 'DMSO']['n_neg'])
mean_n_neg = np.mean(data_screen[data_screen['treatment'] == 'DMSO']['n_neg'])
print(mean_n_neg-3*std_n_neg)
print(mean_n_neg+3*std_n_neg)
data_flt = data_screen[(data_screen['n_pos'] < mean_n_pos-3*std_n_pos) |
                       (data_screen['n_pos'] > mean_n_pos+3*std_n_pos) |
                       (data_screen['n_neg'] < mean_n_neg-3*std_n_neg) |
                       (data_screen['n_neg'] > mean_n_pos+3*std_n_neg) |
                       (data_screen['HSR/DM'] < (mean_ratio-3*std_ratio)) |
                       (data_screen['HSR/DM'] > (mean_ratio+3*std_ratio))].copy().reset_index(drop=True)
print(data_flt['treatment'])
for i in range(len(data_flt)):
    plt.text(x=data_flt['n_neg'][i]+0.005, y=data_flt['n_pos'][i]+0.005, s=data_flt['treatment'][i], size=7, color=(0/255, 191/255, 255/255))
xdata = np.array(data_screen[data_screen['treatment'] == 'DMSO']['n_neg'])
ydata = np.array(data_screen[data_screen['treatment'] == 'DMSO']['n_pos'])
popt, pcov = curve_fit(func, xdata, ydata)
xplot = np.linspace(0, 2500, 100)
plt.plot(xplot, func(xplot, *popt), '--', label='fit: a=%5.3f' % tuple(popt), color=ctrl_colors[0])
x1, y1 = [0, 2500], [0, 2500]
plt.plot(x1, y1, '--', color='r')
plt.xlim([0, 1800])
plt.ylim([0, 1800])
plt.legend()
plt.savefig('%s/R%s_%s_n-neg_vs_n-pos.pdf' % (output_dir, rep, cell))
plt.show()








