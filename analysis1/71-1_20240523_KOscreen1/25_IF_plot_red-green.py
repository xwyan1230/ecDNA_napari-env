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
rep = 2
GFP_cutoff = [3.6, 3.0, 3.2, 2.9]
mCherry_cutoff = [3.1, 3.7, 2.9, 3.4]

data = pd.read_csv('%s/%s/24_%s_IF_rep%s.txt' % (output_dir, sample, sample, rep), na_values=['.'], sep='\t')

fig, ax = plt.subplots(figsize=(9, 6))
fig.subplots_adjust(right=0.8)
sns.histplot(data=data, x='hoechst')
# plt.axvline(hc[0], 0, 1000, c='r')
# plt.axvline(hc[1], 0, 1000, c='r')
# if not os.path.exists("%s%s/%s_%s/hoechst_hist/" % (output_dir, batch, batch, plate)):
#     os.makedirs("%s%s/%s_%s/hoechst_hist/" % (output_dir, batch, batch, plate))
# plt.savefig('%s%s/%s_%s/hoechst_hist/%s.pdf' % (output_dir, batch, batch, plate, sample))
plt.show()

data['log10_GFP'] = np.log10(data['GFP'])
data['log10_mCherry'] = np.log10(data['mCherry'])
data['group'] = ['NA'] * len(data)
data.loc[((data['log10_GFP'] > GFP_cutoff[0]) & (data['log10_mCherry'] < GFP_cutoff[1])) | ((data['log10_GFP'] > GFP_cutoff[2]) & (data['log10_mCherry'] < GFP_cutoff[3])), 'group'] = 'GFP'
data.loc[((data['log10_mCherry'] > mCherry_cutoff[0]) & (data['log10_GFP'] < mCherry_cutoff[1])) | ((data['log10_mCherry'] > mCherry_cutoff[2]) & (data['log10_GFP'] < mCherry_cutoff[3])), 'group'] = 'mCherry'

print(data.head())

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data, x='log10_mCherry', y='log10_GFP', alpha=0.5, s=30)
sns.scatterplot(data=data[data['group'] == 'mCherry'], x='log10_mCherry', y='log10_GFP', color='r', alpha=0.5, s=30)
sns.scatterplot(data=data[data['group'] == 'GFP'], x='log10_mCherry', y='log10_GFP', color='g', alpha=0.5, s=30)
plt.savefig('%s/%s/25_%s_IF_rep%s_red-green.pdf' % (output_dir, sample, sample, rep))
plt.xlim([2.5, 5])
plt.ylim([3.2, 5])
plt.savefig('%s/%s/25_%s_IF_rep%s_red-green_fixscale.pdf' % (output_dir, sample, sample, rep))
plt.show()

print(len(data[data['group'] == 'GFP']))
print(len(data[data['group'] == 'mCherry']))

data.to_csv('%s/%s/25_%s_IF_rep%s_red_green_group.txt' % (output_dir, sample, sample, rep), index=False, sep='\t')

if os.path.exists("%s/IF_red_green.txt" % output_dir):
    pd_group = pd.read_csv('%s/IF_red_green.txt' % output_dir, na_values=['.'], sep='\t')
else:
    pd_group = pd.DataFrame(columns=['sample', 'rep', 'GFP_cutoff_0', 'GFP_cutoff_1', 'GFP_cutoff_2', 'GFP_cutoff_3', 'mCherry_cutoff_0', 'mCherry_cutoff_1', 'mCherry_cutoff_2', 'mCherry_cutoff_3'])
if len(pd_group[(pd_group['sample'] == sample) & (pd_group['rep'] == rep)]) == 0:
    pd_group.loc[len(pd_group.index)] = [sample, rep, GFP_cutoff[0], GFP_cutoff[1], GFP_cutoff[2], GFP_cutoff[3], mCherry_cutoff[0], mCherry_cutoff[1], mCherry_cutoff[2], mCherry_cutoff[3]]
else:
    location = pd_group[(pd_group['sample'] == sample) & (pd_group['rep'] == rep)].index
    pd_group.loc[location] = [sample, rep, GFP_cutoff[0], GFP_cutoff[1], GFP_cutoff[2], GFP_cutoff[3], mCherry_cutoff[0], mCherry_cutoff[1], mCherry_cutoff[2], mCherry_cutoff[3]]
pd_group.to_csv('%s/IF_red_green.txt' % output_dir, index=False, sep='\t')

print("DONE!")