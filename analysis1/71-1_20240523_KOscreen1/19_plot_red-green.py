import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import imutils
import pandas as pd
import shared.image as ima
import tifffile as tif
import math
import nd2
import shared.segmentation as seg
import shared.objects as obj
import os
import seaborn as sns
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'G2'
samples = ['G2']
total_fov = 16
GFP_cutoff = [3.4, 2.9]
mCherry_cutoff = [3, 3.55]

data = pd.DataFrame()
for s in samples:
    df = pd.read_csv('%s/%s/14_%s_red_green.txt' % (output_dir, s, s), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0).reset_index(drop=True)
data['log10_GFP'] = np.log10(data['GFP'])
data['log10_mCherry'] = np.log10(data['mCherry'])
data['group'] = ['NA'] * len(data)
data.loc[(data['log10_GFP'] > GFP_cutoff[0]) & (data['log10_mCherry'] < GFP_cutoff[1]), 'group'] = 'GFP'
data.loc[(data['log10_mCherry'] > mCherry_cutoff[0]) & (data['log10_GFP'] < mCherry_cutoff[1]), 'group'] = 'mCherry'

print(data.head())

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data, x='log10_mCherry', y='log10_GFP', alpha=0.5, s=30)
sns.scatterplot(data=data[data['group'] == 'mCherry'], x='log10_mCherry', y='log10_GFP', color='r', alpha=0.5, s=30)
sns.scatterplot(data=data[data['group'] == 'GFP'], x='log10_mCherry', y='log10_GFP', color='g', alpha=0.5, s=30)
plt.savefig('%s/%s/19_%s_red-green.pdf' % (output_dir, sample, sample))
plt.xlim([2.5, 5])
plt.ylim([3.3, 5])
plt.savefig('%s/%s/19_%s_red-green_fixscale.pdf' % (output_dir, sample, sample))
plt.show()

data.to_csv('%s/%s/19_%s_red_green_group.txt' % (output_dir, sample, sample), index=False, sep='\t')

if os.path.exists("%s/red_green.txt" % output_dir):
    pd_group = pd.read_csv('%s/red_green.txt' % output_dir, na_values=['.'], sep='\t')
else:
    pd_group = pd.DataFrame(columns=['sample', 'GFP_cutoff_0', 'GFP_cutoff_1', 'mCherry_cutoff_0', 'mCherry_cutoff_1'])
if len(pd_group[pd_group['sample'] == sample]) == 0:
    pd_group.loc[len(pd_group.index)] = [sample, GFP_cutoff[0], GFP_cutoff[1], mCherry_cutoff[0], mCherry_cutoff[1]]
else:
    location = pd_group[pd_group['sample'] == sample].index
    pd_group.loc[location] = [sample, GFP_cutoff[0], GFP_cutoff[1], mCherry_cutoff[0], mCherry_cutoff[1]]
pd_group.to_csv('%s/red_green.txt' % output_dir, index=False, sep='\t')

print("DONE!")