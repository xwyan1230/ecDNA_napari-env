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

folders = ['48hr_384well']

# 24hr_384well
"""test_lst = [['XY0%s' % (x+3) for x in range(7)] + ['XY%s' % (x+10) for x in range(3)],
             ['XY%s' % (x+13) for x in range(10)], ['XY%s' % (x+23) for x in range(10)],
             ['XY%s' % (x+33) for x in range(10)], ['XY%s' % (x+43) for x in range(10)],
             ['XY%s' % (x+53) for x in range(10)], ['XY%s' % (x+63) for x in range(10)]]"""
# 48hr_384well
test_lst = [['XY%s' % (x+13) for x in range(10)], ['XY%s' % (x+23) for x in range(10)],
             ['XY%s' % (x+33) for x in range(10)], ['XY%s' % (x+43) for x in range(10)],
             ['XY%s' % (x+53) for x in range(10)], ['XY%s' % (x+63) for x in range(10)]]

ctrl_range = ([0, 500], [500, 1000], [1000, 1500], [1500, 2000], [2000, 2500], [2500, 3000], [3000, 3500], [3500, 4000], [4000, 4500], [4500, 5000])
ctrl_y = [250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750]

data = pd.DataFrame()
for folder in folders:
    df = pd.read_csv('%s/%s/%s_update.txt' % (data_dir, folder, folder), na_values=['.'], sep='\t')
    df['treatment'] = [df['sample_name'][x].split('_')[0] for x in range(len(df))]
    df['HSR/DM'] = (df['n_pos']+0.01)/(df['n_neg']+0.01)
    data = pd.concat([data, df], axis=0)

data_ctrl = data[data['treatment'] == 'DMSO'].copy().reset_index(drop=True)
ctrl_x = []
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
plt.savefig('%s/%s_ctrl_mean.pdf' % (output_dir, folders[0]))
plt.show()