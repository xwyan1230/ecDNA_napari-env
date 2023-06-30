import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.dataframe as dat
import shared.math as mat
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230504_analysis_DM_6hr/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample_name = 'G8'
combine_lst = ['G8_1', 'G8_2']

data = pd.DataFrame()
for sample in combine_lst:
    df = pd.read_csv('%s/txt/%s_n4.txt' % (data_dir1, sample), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0)

data.to_csv('%s/txt/%s_n4.txt' % (output_dir, sample_name), index=False, sep='\t')
print("DONE!")