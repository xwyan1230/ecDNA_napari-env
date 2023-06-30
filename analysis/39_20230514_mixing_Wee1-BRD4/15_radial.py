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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230514_analysis_mixing_Wee1-BRD4/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

exp = '20230512_mixing_Wee1-BRD4_6hr'
sample = '10_BRD4-GFP-6hr_Ctrl-mCh-6hr'
GFP_sample = 'BRD4_6hr'
mCherry_sample = 'Ctrl'
hue_order = [mCherry_sample, GFP_sample]

df = pd.read_csv('%s%s/%s_summary.txt' % (data_dir1, sample, sample), na_values=['.'], sep='\t')
sample_lst = []
for i in range(len(df)):
    if df['group'][i] == 'GFP':
        sample_lst.append(GFP_sample)
    elif df['group'][i] == 'mCherry':
        sample_lst.append(mCherry_sample)
    else:
        sample_lst.append('NA')
df['sample'] = sample_lst

feature_lst = ['nuclear_int', 'DNAFISH_int', 'DNAFISH_seg_label', 'int_r_to_edge', 'int_r_to_center', 'int_relative_r']
for f in feature_lst:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]
hue_order1 = ['background', 'DNAFISH']


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask) if e != 0]
    out = list(map(img.__getitem__, index))
    return out


df_r = pd.DataFrame()
# new approach
for sample in hue_order:
    print(sample)
    data = df[df['sample'] == sample].copy().reset_index(drop=True)
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
df_r = pd.concat([df_r, df_r_temp], axis=0)
if not os.path.exists("%s%s/txt_radial/" % (output_dir, sample)):
    os.makedirs("%s%s/txt_radial/" % (output_dir, sample))
df_r.to_csv('%s%s/txt_radial/radial.txt' % (output_dir, sample), index=False, sep='\t')
print("DONE!")
