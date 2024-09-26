import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
from shared.sinaplot import sinaplot
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.dataframe as dat
import shared.math as mat
import math
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import iqr
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240222_analysis_DF-HFH_density/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

folders = ['24hr_384well']
test_lst = [['XY0%s' % (x+3) for x in range(7)] + ['XY%s' % (x+10) for x in range(3)],
             ['XY%s' % (x+13) for x in range(10)], ['XY%s' % (x+23) for x in range(10)],
             ['XY%s' % (x+33) for x in range(10)], ['XY%s' % (x+43) for x in range(10)],
             ['XY%s' % (x+53) for x in range(10)], ['XY%s' % (x+63) for x in range(10)]]

sample_lst = ['XY0%s' % (x+3) for x in range(7)] + ['XY%s' % (x+10) for x in range(63)]  # 24hr_384well
# sample_lst = ['XY%s' % (x+11) for x in range(62)]  # 48hr_384well
# sample_lst = ['XY11', 'XY12', 'XY13', 'XY14'] + ['XY%s' % (x+31) for x in range(42)]  # 72hr_384well

density_lst = [10000, 8000, 6000, 4000, 3000, 2000, 1000]
# size_lst = [40, 35,30,25,20,15,10]
ctrl_colors = [(0.6, 0.81, 0.98), (0.5, 0.81, 0.98), (0.4, 0.81, 0.98), (0.3, 0.81, 0.98), (0.2, 0.81, 0.98), (0.1, 0.81, 0.98), (0, 0.81, 0.98)][::-1]
MK1775_colors = [(0.88, 0.6, 0.6), (0.88, 0.5, 0.5), (0.88, 0.4, 0.4), (0.88, 0.3, 0.3), (0.88, 0.2, 0.2), (0.88, 0.1, 0.1), (0.88, 0, 0)][::-1]
rainboo_colors = [(220/255, 20/255, 60/255), (1, 140/255, 0), (1, 215/255, 0), (154/255, 205/255, 50/255), (64/255, 224/255, 208/255),
                  (100/255, 149/255, 237/255), (123/255, 104/255, 238/255)]

"""ctrl_range = ([500, 1000], [1000, 1500], [1500, 2000], [2000, 2500], [2500, 3000], [3000, 4000])
ctrl_y = [750, 1250, 1750, 2250, 2750, 3500]"""

data = pd.DataFrame()
for folder in folders:
    df = pd.read_csv('%s/%s/%s_update1_hc.txt' % (data_dir, folder, folder), na_values=['.'], sep='\t')
    df['treatment'] = [df['sample_name'][x].split('_')[0] for x in range(len(df))]
    df['HSR/DM'] = (df['n_pos']+0.01)/(df['n_neg']+0.01)
    data = pd.concat([data, df], axis=0)

data = data[data['sample'].isin(sample_lst)].copy().reset_index(drop=True)
data_ctrl = data[data['treatment'] == 'DMSO'].copy().reset_index(drop=True)
"""ctrl_x = []
for i in ctrl_range:
    temp = data_ctrl[(data_ctrl['n_total'] >= i[0]) & (data_ctrl['n_total'] < i[1])]
    if len(temp) > 0:
        ctrl_x.append(np.mean(temp['HSR/DM']))
    else:
        ctrl_x.append(-1)

plt.subplots(figsize=(9, 6))
feature = 'HSR/DM'
sns.scatterplot(data=data_ctrl, y=feature, x='n_total', color=(0.4, 0.81, 0.98), alpha=0.8, s=20)
sns.scatterplot(y=ctrl_x, x=ctrl_y, color='r', s=30)
plt.xlim([0, 5000])
plt.ylim([0, 3])
plt.savefig('%s/%s_ctrl_mean_fit.pdf' % (output_dir, folders[0]))
plt.show()"""


def func(x, a, b):
    return a * x + b * x**2 + 1


data_fit = data_ctrl[data_ctrl['n_total'] > 500].copy().reset_index(drop=True)
xdata = np.array(data_fit['n_total'])
ydata = np.array(data_fit['HSR/DM'])
popt, pcov = curve_fit(func, xdata, ydata)
print(popt)
xplot = np.linspace(0, 5000, 100)

plt.subplots(figsize=(9, 6))
feature = 'HSR/DM'
sns.scatterplot(data=data_ctrl, y=feature, x='n_total', color=(0.4, 0.81, 0.98), alpha=0.8, s=20)
# sns.scatterplot(y=ctrl_x, x=ctrl_y, color='r', s=30)
plt.plot(xplot, func(xplot, *popt), '--', label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
plt.legend()
plt.xlim([0, 5000])
plt.ylim([0, 3])
plt.savefig('%s/%s_ctrl_mean_fit_curve_update1_hc.pdf' % (output_dir, folders[0]))
plt.show()

data['log2FC_HSR-DM'] = [np.log2(data['HSR/DM'][j]/func(data['n_total'][j], *popt)) for j in range(len(data))]

plt.subplots(figsize=(9, 6))
feature = 'log2FC_HSR-DM'
for i in range(len(test_lst)):
    test = test_lst[i]
    data_screen = data[data['sample'].isin(test)].copy().reset_index(drop=True)

    sns.scatterplot(data=data_screen, x=feature, y='n_total', alpha=0.8, s=30, color=MK1775_colors[i])
    sns.scatterplot(data=data_screen[data_screen['treatment'] == 'DMSO'], x=feature, y='n_total', color=ctrl_colors[i], alpha=0.8, s=30)
plt.xlim([-1.5, 1.5])
plt.ylim([0, 5000])
plt.savefig('%s/%s_%s_MK1775_update_log2FC-corrected_update1_hc.pdf' % (output_dir, folder, feature))
plt.show()

plt.subplots(figsize=(9, 6))
feature = 'log2FC_HSR-DM'
for i in range(len(test_lst)):
    test = test_lst[i]
    data_screen = data[data['sample'].isin(test)].copy().reset_index(drop=True)

    sns.scatterplot(data=data_screen, x=feature, y='n_total', alpha=0.8, s=30, color=rainboo_colors[i])
plt.xlim([-1.5, 1.5])
plt.ylim([0, 5000])
plt.savefig('%s/%s_%s_density_update_log2FC-corrected_update1_hc.pdf' % (output_dir, folder, feature))
plt.show()







