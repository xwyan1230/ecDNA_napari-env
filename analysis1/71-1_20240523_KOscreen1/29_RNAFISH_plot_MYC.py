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
from skimage.measure import label, regionprops
import seaborn as sns
from scipy.stats import ttest_ind

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'G2'
rep = 2
hueorder = ['mCherry', 'GFP']

data = pd.read_csv('%s/%s/28_%s_RNAFISH_rep%s_red_green_group.txt' % (output_dir, sample, sample, rep), na_values=['.'], sep='\t')
data['log10_RNAFISH'] = np.log10(data['RNAFISH'])
data_flt = data.sort_values(by='group', ascending=False).copy().reset_index(drop=True)

fig, ax = plt.subplots(figsize=(4, 6))
fig.subplots_adjust(left=0.4)
ax = sns.boxplot(data=data_flt, x='group', y='log10_RNAFISH', hue_order=hueorder)
lines = ax.get_lines()
categories = ax.get_xticks()
plt.close()

fig, ax = plt.subplots(figsize=(4, 6))
fig.subplots_adjust(left=0.4)
sns.swarmplot(data=data_flt, x='group', y='log10_RNAFISH', hue_order=hueorder, s=1, color=(255 / 255, 140 / 255, 0 / 255))

nobs = [len(data_flt[data_flt['group'] == hueorder[i]]) for i in range(len(hueorder))]
nobs = ["%s" % str(x) for x in nobs]
for cat in categories:
    # every 4th line at the interval of 6 is median line
    # 0 -> p25 1 -> p75 2 -> lower whisker 3 -> upper whisker 4 -> p50 5 -> upper extreme value
    y = round(lines[4 + cat * 6].get_ydata()[0], 2)
    ax.text(cat, y, '%s (%s)' % (y, nobs[cat]), ha='center', va='center', fontweight='bold', size=7, color='white', bbox=dict(facecolor='#445A64'))

p = ttest_ind(data_flt[data_flt['group'] == 'mCherry']['log10_RNAFISH'].tolist(), data_flt[data_flt['group'] == 'GFP']['log10_RNAFISH'].tolist())
print(p)
"""y, h, col = data_flt['log10_DNAFISH_total_int_merge'].max() + 0.05, 0.05, 'k'
plt.plot([0, 0, 1, 1], [y, y + h, y + h, y], lw=1, c=col)
plt.text((0 + 1) * .5, y + h, "ns", ha='center', va='bottom', color=col)"""
plt.savefig('%s/%s/29_%s_RNAFISH_rep%s_MYC.pdf' % (output_dir, sample, sample, rep))
plt.ylim([2.7, 5])
plt.savefig('%s/%s/29_%s_RNAFISH_rep%s_MYC_fixscale.pdf' % (output_dir, sample, sample, rep))
plt.show()

print("DONE!")
