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

exp = '20230409_DM_plate_DNAFISH_control'

feature_lst = ['nuclear_int', 'DNAFISH_int', 'DNAFISH_seg_label', 'int_r_to_edge', 'int_r_to_center', 'int_relative_r']
hue_order1 = ['background', 'DNAFISH']
row_lst = ['B', 'C', 'D', 'E', 'F', 'G']
column_lst = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
"""sample_lst = []
for row in row_lst:
    sample_lst = sample_lst + ['%s%s' % (row, column) for column in column_lst]"""


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask) if e != 0]
    out = list(map(img.__getitem__, index))
    return out


# new approach
for r in range(len(row_lst)):
    for c in range(len(column_lst)):
        sample = '%s%s' % (row_lst[r], column_lst[c])
        print(sample)
        if os.path.exists('%s/txt/%s_n4.txt' % (data_dir1, sample)):
            df = pd.read_csv('%s/txt/%s_n4.txt' % (data_dir1, sample), na_values=['.'], sep='\t')
            for f in feature_lst:
                df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]
            data = df
            df_r_temp = pd.DataFrame()
            sample_category_lst = []
            intensity_lst = []
            seg_lst = []
            relative_r_lst = []
            for i in range(len(data)):
                sample_category_lst = sample_category_lst + [hue_order1[0]] * len(data['int_relative_r'][i]) + [
                    hue_order1[1]] * len(data['int_relative_r'][i])
                relative_r_lst = relative_r_lst + data['int_relative_r'][i] + data['int_relative_r'][i]
                intensity_lst = intensity_lst + data['nuclear_int'][i] + data['DNAFISH_int'][i]
                seg_lst = seg_lst + [1] * len(data['int_relative_r'][i]) + data['DNAFISH_seg_label'][i]
            df_r_temp['sample_category'] = sample_category_lst
            df_r_temp['relative_r'] = relative_r_lst
            df_r_temp['intensity'] = intensity_lst
            df_r_temp['seg'] = seg_lst
            df_r_temp['sample'] = [sample] * len(df_r_temp)
            len_bg = len(df_r_temp[df_r_temp['sample_category'] == hue_order1[0]])
            len_sample = len(df_r_temp[(df_r_temp['sample_category'] == hue_order1[1]) & (df_r_temp['seg'] == 1)])
            df_r_temp['weights'] = [1.0 / len_bg if i == hue_order1[0] else 1.0 / len_sample for i in
                                    df_r_temp['sample_category']]
            if not os.path.exists("%s/txt_radial/" % output_dir):
                os.makedirs("%s/txt_radial/" % output_dir)
            df_r_temp.to_csv('%s/txt_radial/%s_radial.txt' % (output_dir, sample), index=False, sep='\t')
