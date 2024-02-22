import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import shared.dataframe as dat
import shared.objects as obj
import shared.math as mat
import cv2
import math
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230724_analysis_Natasha_colcemid/TPL/mCh-DMSO_GFP-TPL/"
output_dir = "%s06_analysis/" % master_folder

sample = 'mCh-DMSO_GFP-TPL'

df = pd.read_csv('%s/06_analysis/%s_3_n2.txt' % (master_folder, sample), na_values=['.'], sep='\t')
df['ln_mean_int_green'] = np.log(df['mean_int_green'])
df['ln_mean_int_red'] = np.log(df['mean_int_red'])

plt.subplots(figsize=(6, 6))
sns.scatterplot(data=df, y='ln_mean_int_green', x='ln_mean_int_red', color=(0.30, 0.30, 0.30), s=5, alpha=0.5)
# plt.savefig('%s/%s/%s_red-green_scatter.pdf' % (output_dir, sample, sample))
plt.show()