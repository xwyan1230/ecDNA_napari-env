import skimage.io as skio
import napari
import imutils
import shared.image as ima
import shared.objects as obj
from skimage.morphology import disk, dilation, medial_axis
import shared.dataframe as dat
import math
import shared.display as dis
from skimage.morphology import disk, dilation
import tifffile as tif
import matplotlib.pyplot as plt
import seaborn as sns
from skimage.measure import label, regionprops
import numpy as np
import pandas as pd
# from napari_animation import Animation
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230706_analysis_radial/"
data_dir = "%sdata/raw/" % master_folder
data_dir1 = "%sdata/seg_tif/" % master_folder
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

# sample = 'Colo320DM_acidFISH_lamin_3d'
sample = 'Colo320HSR_acidFISH_lamin_3d'
# prefix = '20230614_acidFISH_lamin_ColoDM_DM'
prefix = '20230614_acidFISH_lamin_ColoHSR_HSR'

df = pd.read_csv('%s%s/%s_n4_full.txt' % (data_dir2, sample, sample), na_values=['.'], sep='\t')
n_nuclear_convex_dilation = 4
cn_lst = []
for i in range(len(df)):
    x = np.log(df['total_int_DNAFISH_proj'][i])
    if x < 17:
        cn_lst.append('<17')
    elif x < 17.5:
        cn_lst.append('[17, 17.5)')
    elif x < 18:
        cn_lst.append('[17.5, 18)')
    elif x < 18.5:
        cn_lst.append('[18, 18.5)')
    elif x >= 18.5:
        cn_lst.append('>18.5')
    else:
        cn_lst.append('NA')
        print(df['total_int_DNAFISH_proj'][i])
        print("NA")

df['cn_group'] = cn_lst
df.to_csv('%s%s/%s_n%s_full.txt' % (output_dir, sample, sample, n_nuclear_convex_dilation), index=False, sep='\t')

print("DONE!")