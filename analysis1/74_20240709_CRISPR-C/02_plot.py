import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
from skimage.measure import regionprops, label
from skimage.morphology import binary_dilation, disk
import numpy as np
import shared.segmentation as seg
import shared.objects as obj
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240709_analysis_CRISPR-C/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

hc = [5000, 40000]
cutoff = 3.65

sample = 'XY06'

df = pd.DataFrame(columns=['sample', 'identity', 'n_filtered', 'percentage'])
data = pd.read_csv('%s/txt/%s.txt' % (data_dir, sample), na_values=['.'], sep='\t')
data['log10_GFP'] = np.log10(data['GFP'])

plt.subplots(figsize=(9, 7))
plt.hist(data['hoechst'], weights=np.ones(len(data)) / len(data), color='w', edgecolor='black')
plt.show()

data_filtered = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])]
n_filtered = len(data_filtered)
n_pos = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_GFP'] > cutoff)])
percent = n_pos/n_filtered
print(percent)

plt.subplots(figsize=(9, 7))
plt.hist(data_filtered['log10_GFP'], weights=np.ones(len(data_filtered)) / len(data_filtered), range=[2.5, 5], bins=40, color='w', edgecolor='black')
plt.axvline(x=cutoff, color='#bc4d4a', linestyle='--')
# plt.xlim([2.5, 5])
# plt.ylim([0, 0.6])
# if not os.path.exists('%s/%s/' % (output_dir, sample)):
#     os.makedirs('%s/%s/' % (output_dir, sample))
# plt.savefig('%s/%s/%s_cutoff_qc_%s.pdf' % (output_dir, sample, sample))
plt.show()
