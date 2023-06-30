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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

hue_order = ['G1', 'S', 'G2', 'M']
line_colors = [(220/255, 20/255, 60/255), (154/255, 205/255, 50/255), (255/255, 127/255, 80/255), (0.85, 0.35, 0.25)]

sample = 'DM_3_49pos'
# sample = 'DM_324pos_merge'
df = pd.read_csv('%s%s_summary.txt' % (output_dir, sample), na_values=['.'], sep='\t')
# df = pd.read_csv('%s%s_cellcycle_updated1.txt' % (output_dir, sample), na_values=['.'], sep='\t')
df['total_int_DNAFISH'] = df['area_nuclear']*df['mean_int_DNAFISH']
df['total_int_hoechst'] = df['area_nuclear_IF']*df['mean_int_hoechst']
df['total_int_ecDNA'] = df['total_area_ecDNA']*df['mean_int_ecDNA']
df['ln_mean_int_green'] = np.log(df['mean_int_green'])
df['ln_mean_int_red'] = np.log(df['mean_int_red'])
df['ln_mean_int_red_cal'] = df['ln_mean_int_red'] - 0.14 * df['ln_mean_int_green']
df['ln_mean_int_green_cal'] = df['ln_mean_int_green'] - 0.2 * df['ln_mean_int_red']
# df['total_int_hoechst'] = df['area_nuclear']*df['mean_int_nuclear']

"""sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
feature = 'mean_int_DNAFISH'
sinaplot(data=df, x='cellcycle', y=feature, order=hue_order, violin=False, scale='area', point_size=2)
plt.ylim([0, 20000])
plt.savefig('%s/%s_%s.pdf' % (output_dir, sample, feature))
plt.show()"""

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='ln_mean_int_green', x='ln_mean_int_red', c=df['total_int_ecDNA'], cmap='Reds', vmin=0, s=10)
# sns.scatterplot(data=df, y='ln_mean_int_green', x='ln_mean_int_red', color=(0.30, 0.30, 0.30), s=5, alpha=0.5)
# sns.scatterplot(data=df[df['cellcycle'] == 'G2'], y='ln_mean_int_green', x='ln_mean_int_red', color='red', s=5)
plt.show()