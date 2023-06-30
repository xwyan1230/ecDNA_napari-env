import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

n_dilation = 4

# samples
sample = 'DM_3_49pos'
figure_name = 'DM'
df = pd.read_csv('%s%s_summary.txt' % (data_dir, sample), na_values=['.'], sep='\t')

hue_order_r = ['17.375', '17.625', '17.875', '18.125', '18.375', '18.625', '18.875']
hue_order1 = ['background', 'DNAFISH']

feature = ['nuclear_int', 'DNAFISH_int', 'DNAFISH_seg_label', 'int_r_to_edge', 'int_r_to_center', 'int_relative_r']

for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

df_sample = df[df['cellcycle'].isin(['G1', 'S', 'G2'])].copy().reset_index(drop=True)
df = df_sample
df['total_int_MYC'] = df['area_nuclear_IF'] * df['mean_int_MYC']
df['ln_total_int_MYC'] = np.log(df['total_int_MYC'])
print(len(df))


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask) if e != 0]
    out = list(map(img.__getitem__, index))
    return out


# new approach
df_r = pd.DataFrame()
for k in range(len(hue_order_r)):
    cutoff = float(hue_order_r[k])
    data = df[(df['ln_total_int_MYC'] >= cutoff - 0.125) & (df['ln_total_int_MYC'] < cutoff + 0.125)].copy().reset_index(drop=True)
    df_r_temp = pd.DataFrame()
    sample_category_lst = []
    intensity_lst = []
    seg_lst = []
    relative_r_lst = []
    for i in range(len(data)):
        sample_category_lst = sample_category_lst + [hue_order1[0]] * len(data['int_relative_r'][i]) + [hue_order1[1]] * len(data['int_relative_r'][i])
        relative_r_lst = relative_r_lst + data['int_relative_r'][i] + data['int_relative_r'][i]
        intensity_lst = intensity_lst + data['nuclear_int'][i] + data['DNAFISH_int'][i]
        seg_lst = seg_lst + [1] * len(data['int_relative_r'][i]) + data['DNAFISH_seg_label'][i]
    df_r_temp['sample_category'] = sample_category_lst
    df_r_temp['relative_r'] = relative_r_lst
    df_r_temp['intensity'] = intensity_lst
    df_r_temp['seg'] = seg_lst
    df_r_temp['sample'] = [hue_order_r[k]] * len(df_r_temp)
    len_bg = len(df_r_temp[df_r_temp['sample_category'] == hue_order1[0]])
    len_sample = len(df_r_temp[(df_r_temp['sample_category'] == hue_order1[1]) & (df_r_temp['seg'] == 1)])
    df_r_temp['weights'] = [1.0/len_bg if i == hue_order1[0] else 1.0/len_sample for i in df_r_temp['sample_category']]
    df_r = pd.concat([df_r, df_r_temp], axis=0)

df_r.to_csv('%s%s_radial_MYC.txt' % (output_dir, sample), index=False, sep='\t')

print("DONE!")