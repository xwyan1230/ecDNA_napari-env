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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230706_analysis_radial/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

# sample_name = 'Colo320DM_acidFISH_lamin_3d'
sample_name = 'Colo320HSR_acidFISH_lamin_3d'
hue_order = ['<17', '[17, 17.5)', '[17.5, 18)', '[18, 18.5)', '>18.5']

df = pd.read_csv('%s%s/%s_n4_full.txt' % (data_dir, sample_name, sample_name), na_values=['.'], sep='\t')

feature_lst = ['DNAFISH_seg_label', 'int_r_to_edge', 'int_relative_r']
for f in feature_lst:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

df_r = pd.DataFrame()
# new approach
for sample in hue_order:
    print(sample)
    data = df[df['cn_group'] == sample].copy().reset_index(drop=True)
    df_r_temp = pd.DataFrame()
    seg_lst = []
    relative_r_lst = []
    absolute_r_lst = []
    for i in range(len(data)):
        relative_r_lst = relative_r_lst + data['int_relative_r'][i]
        seg_lst = seg_lst + data['DNAFISH_seg_label'][i]
        absolute_r_lst = absolute_r_lst + data['int_r_to_edge'][i]
    df_r_temp['relative_r'] = relative_r_lst
    df_r_temp['seg'] = seg_lst
    df_r_temp['absolute_r'] = absolute_r_lst
    df_r_temp['sample'] = [sample] * len(df_r_temp)
    len_sample = len(df_r_temp[df_r_temp['seg'] == 1.0])
    print(len_sample)
    # df_r_temp['weights'] = [1.0 / len_sample] * len(df_r_temp)
    df_r = pd.concat([df_r, df_r_temp], axis=0)
if not os.path.exists("%s%s/" % (output_dir, sample_name)):
    os.makedirs("%s%s/" % (output_dir, sample_name))
df_r.to_csv('%s%s/%s_radial_CN.txt' % (output_dir, sample_name, sample_name), index=False, sep='\t')
print("DONE!")
