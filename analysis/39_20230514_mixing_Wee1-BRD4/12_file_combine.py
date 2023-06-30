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
import shared.math as mat
import math
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230514_analysis_mixing_Wee1-BRD4/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

exp = '20230512_mixing_Wee1-BRD4_6hr'
sample = '9_Ctrl-GFP-6hr_BRD4-mCh-6hr'
"""file_name = 'alignment_label'
file1 = 'alignment_label_1'
file2 = 'alignment_label_2'"""
file_name = '%s_n4' % sample
file1 = '%s_1_n4' % sample
file2 = '%s_2_n4' % sample
df1 = pd.read_csv('%s%s/%s.txt' % (data_dir1, sample, file1), na_values=['.'], sep='\t')
df2 = pd.read_csv('%s%s/%s.txt' % (data_dir1, sample, file2), na_values=['.'], sep='\t')
df = pd.concat([df1, df2], axis=0)
df.to_csv('%s%s/%s.txt' % (output_dir, sample, file_name), index=False, sep='\t')
print("DONE!")