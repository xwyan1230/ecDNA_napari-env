import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
# from shared.sinaplot import sinaplot
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
from matplotlib_venn import venn2
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

folder = '5uM_48hr'
batches = [1, 2, 3]
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
colors = ['#bc4d4a', '#669daa', '#dd933c'] * 3
fov_area = 2.336  # mm2
well_area = 11
ylim_val = [-500, 17500]

data = pd.DataFrame()
for batch in batches:
    df = pd.read_csv('%s/%s/%s_%s_update.txt' % (data_dir, folder, folder, batch), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0)

data = data[data['n_filtered']>500].copy().reset_index(drop=True)

"""plt.subplots(figsize=(9, 7))
m=sns.histplot(data=data['log2_h2-3'], bins=40)
plt.show()"""

"""plt.subplots(figsize=(9, 7))
sns.histplot(data=data['green_median'], bins=40)
plt.show()

plt.subplots(figsize=(9, 7))
sns.histplot(data=data['red_median'], bins=40)
plt.show()

plt.subplots(figsize=(9, 7))
sns.histplot(data=data['per_green'], bins=40)
plt.show()

plt.subplots(figsize=(9, 7))
sns.histplot(data=data['per_red'], bins=40)
plt.show()"""

print(data[((data['per_green']<0.1) & (data['per_red']<0.1))|(data['per_NA']>0.5)]['plate'])
print(data[((data['per_green']<0.1) & (data['per_red']<0.1))|(data['per_NA']>0.5)]['sample'])
print(data[((data['per_green']<0.1) & (data['per_red']<0.1))|(data['per_NA']>0.5)]['per_NA'])
# print(data[(data['per_green']<0.15) & (data['per_red']<0.15)]['red_median'])

"""23     1
156    1
313    2
359    2
366    2
443    3
Name: plate, dtype: int64
23      XY24
156    XY165
313    XY137
359    XY185
366    XY192
443     XY76
Name: sample, dtype: object
23     12.123327
156    11.838022
313    12.171077
359    12.748049
366    12.880945
443    11.835092
Name: green_median, dtype: float64
23     13.248995
156    14.018687
313    13.575263
359    13.531347
366    12.197406
443    14.069655
Name: red_median, dtype: float64"""








