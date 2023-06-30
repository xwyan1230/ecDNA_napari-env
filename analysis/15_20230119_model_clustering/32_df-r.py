import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import matplotlib as mpl
import numpy as np
import math
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230119_model_clustering/"
data_dir = "%stxt/dataset6_dataset3-radial/" % master_folder
output_dir = "%stxt/dataset6_dataset3-radial/" % master_folder

samples = ['40-0_5', '40-1_5', '40-2_5', '40-5_5', '40-10_5', '40-25_5', '40-50_5', '40-75_5', '40-100_5', '40-200_5', '40-300_5', '40-400_5', '40-500_5', '40-1000_5',
           '40-2000_5', '40-3000_5', '40-4000_5', '40-5000_5']
# samples = ['0-0_5', '10-0_5', '20-0_5', '30-0_5', '40-0_5', '50-0_5', '60-0_5', '70-0_5']
df_r_name = 'cen-r_40'

cmap = mpl.cm.get_cmap('Spectral')
x = np.arange(0, 1, 1/len(samples))
line_color = [cmap(i) for i in x]
sns.set_palette(sns.color_palette(line_color))

data = pd.DataFrame()
hue_order = samples
for i in samples:
    cen_r = i.split('-')[0]
    sample = i.split('-')[1]
    coefficient = sample.split('_')[0]
    df = pd.read_csv(("%s%s/%s_radial.txt" % (data_dir, cen_r, sample)), na_values=['.'], sep='\t')
    feature = ['DNAFISH_seg_label', 'int_relative_r']
    for f in feature:
        df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]
    data = pd.concat([data, df], axis=0)
df = data

hue_order1 = ['background', 'DNAFISH']

# new approach
df_r = pd.DataFrame()
for k in range(len(hue_order)):
    cen_r = int(hue_order[k].split('-')[0])
    sample = hue_order[k].split('-')[1]
    coefficient = int(sample.split('_')[0])
    print(hue_order[k])
    data = df[(df['cen_r'] == cen_r) & (df['coefficient'] == coefficient)].copy().reset_index(drop=True)
    print(len(data))
    df_r_temp = pd.DataFrame()
    sample_category_lst = []
    seg_lst = []
    relative_r_lst = []
    for i in range(len(data)):
        sample_category_lst = sample_category_lst + [hue_order1[0]] * len(data['int_relative_r'][i]) + [hue_order1[1]] * len(data['int_relative_r'][i])
        relative_r_lst = relative_r_lst + data['int_relative_r'][i] + data['int_relative_r'][i]
        seg_lst = seg_lst + [1] * len(data['int_relative_r'][i]) + data['DNAFISH_seg_label'][i]
    df_r_temp['sample_category'] = sample_category_lst
    df_r_temp['relative_r'] = relative_r_lst
    df_r_temp['seg'] = seg_lst
    df_r_temp['sample'] = [hue_order[k]] * len(df_r_temp)
    len_bg = len(df_r_temp[df_r_temp['sample_category'] == hue_order1[0]])
    len_sample = len(df_r_temp[(df_r_temp['sample_category'] == hue_order1[1]) & (df_r_temp['seg'] == 1)])
    df_r_temp['weights'] = [1.0/len_bg if i == hue_order1[0] else 1.0/len_sample for i in df_r_temp['sample_category']]
    df_r = pd.concat([df_r, df_r_temp], axis=0)
    print(len(df_r))

df_r.to_csv('%sdf_r_%s.txt' % (output_dir, df_r_name), index=False, sep='\t')

print("df_r done!")

