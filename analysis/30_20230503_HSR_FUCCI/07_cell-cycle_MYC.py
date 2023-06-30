import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
import shared.dataframe as dat
from skimage import segmentation
from shared.sinaplot import sinaplot
import shared.math as mat
import math
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230517_analysis_HSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

hue_order = ['G1', 'S', 'G2M']
line_colors = [(220/255, 20/255, 60/255), (154/255, 205/255, 50/255), (255/255, 127/255, 80/255)]

sample = 'HSR_81pos_merge'
df = pd.read_csv('%s%s_cellcycle_hoechstonly.txt' % (output_dir, sample), na_values=['.'], sep='\t')
df['total_int_MYC'] = df['area_nuclear']*df['mean_int_MYC']
df_cc = df[df['cellcycle'].isin(hue_order)]

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
feature = 'total_int_MYC'
sinaplot(data=df, x='cellcycle', y=feature, order=hue_order, violin=False, scale='area', point_size=2)
plt.ylim([0, 3E8])
plt.savefig('%s/%s_%s.pdf' % (output_dir, sample, feature))
plt.show()

# 3123
# 3435
# 2300
"""print(len(df_cc))
print(len(df[df['cellcycle'] == 'G1']))
print(len(df[df['cellcycle'] == 'S']))
print(len(df[df['cellcycle'] == 'G2']))
print(len(df[df['cellcycle'] == 'M']))"""