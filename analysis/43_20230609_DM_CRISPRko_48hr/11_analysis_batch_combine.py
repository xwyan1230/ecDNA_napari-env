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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'C3'
total_batch = 2

df = pd.DataFrame()
for i in range(total_batch):
    df_temp = pd.read_csv('%s/figures/%s/%s_n4_full_%s.txt' % (master_folder, sample, sample, i+1), na_values=['.'], sep='\t')
    df_temp['batch'] = [i+1] * len(df_temp)
    df = pd.concat([df, df_temp], axis=0)

df.to_csv('%s%s/%s_n4_full.txt' % (output_dir, sample, sample), index=False, sep='\t')
print("DONE!")