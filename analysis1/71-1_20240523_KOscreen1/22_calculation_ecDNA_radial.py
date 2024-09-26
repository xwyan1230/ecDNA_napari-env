import skimage.io as skio
import napari
import numpy as np
import pandas as pd
import shared.image as ima
import nd2
import shared.dataframe as dat
from skimage.measure import label, regionprops
from skimage.morphology import medial_axis

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'G2'
samples = ['G2']
dshape_factor = 0.145
pixel_size = 300 / 2720  # uM
local_size = 200
hueorder =['mCherry', 'GFP']

data_rg = pd.read_csv('%s/%s/19_%s_red_green_group.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
data = pd.DataFrame()
for s in samples:
    df = pd.read_csv('%s/%s/18_%s_radial.txt' % (output_dir, s, s), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0).reset_index(drop=True)

data = data.sort_values(['sample', 'fov', 'label_nuclear'], ascending=[True, True, True]).copy().reset_index(drop=True)
data_rg = data_rg.sort_values(['sample', 'fov', 'label_mean_int'], ascending=[True, True, True]).copy().reset_index(drop=True)

if (data['label_nuclear'].tolist() == data_rg['label_mean_int'].tolist()) & (data['fov'].tolist() == data_rg['fov'].tolist()):
    data['group'] = data_rg['group']
    data_flt = data[data['group'].isin(['GFP', 'mCherry'])].copy().reset_index(drop=True)
    feature_lst = ['DNAFISH_seg_label', 'int_r_to_edge', 'int_relative_r']
    for f in feature_lst:
        data_flt[f] = [dat.str_to_float(data_flt[f][i]) for i in range(len(data_flt))]

    df_r = pd.DataFrame()
    for h in range(len(hueorder)):
        df = data_flt[data_flt['group'] == hueorder[h]].copy().reset_index(drop=True)
        # new approach
        df_r_temp = pd.DataFrame()
        seg_lst = []
        relative_r_lst = []
        absolute_r_lst = []
        for i in range(len(df)):
            relative_r_lst = relative_r_lst + df['int_relative_r'][i]
            seg_lst = seg_lst + df['DNAFISH_seg_label'][i]
            absolute_r_lst = absolute_r_lst + df['int_r_to_edge'][i]
        df_r_temp['relative_r'] = relative_r_lst
        df_r_temp['seg'] = seg_lst
        df_r_temp['absolute_r'] = absolute_r_lst
        df_r_temp['group'] = [hueorder[h]] * len(df_r_temp)
        len_sample = len(df_r_temp[df_r_temp['seg'] == 1.0])
        print(len_sample)
        # df_r_temp['weights'] = [1.0 / len_sample] * len(df_r_temp)
        df_r = pd.concat([df_r, df_r_temp], axis=0)
    df_r.to_csv('%s/%s/22_%s_radial_calculated.txt' % (output_dir, sample, sample), index=False, sep='\t')
    
print("DONE!")

