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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230504_analysis_DM_2hr/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

feature_lst = ['sample', 'n', 'averageD_mean', 'averageD_std', 'averageD_0.5',
               'area_nuclear', 'mean_int_nuclear', 'mean_int_DNAFISH', 'n_ecDNA',
               'mean_int_ecDNA', 'total_area_ecDNA', 'total_area_ratio_ecDNA']
row_lst = ['B', 'C', 'D', 'E', 'F', 'G']
column_lst = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
sample_lst = []
for row in row_lst:
    sample_lst = sample_lst + ['%s%s' % (row, column) for column in column_lst]
df_mean = pd.DataFrame(columns=feature_lst)
for sample in sample_lst:
    if os.path.exists('%s/txt/%s_n4.txt' % (data_dir1, sample)):
        df = pd.read_csv('%s/txt/%s_n4.txt' % (data_dir1, sample), na_values=['.'], sep='\t')
        df['r'] = np.sqrt(df['area_nuclear'] / math.pi)
        df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA'] / math.pi)
        df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt'] / df['r']
        df['dis_to_hub_area_normalized'] = df['dis_to_hub_area'] / df['r']
        df_mean.loc[len(df_mean.index)] = [sample, len(df), np.mean(df['dis_to_hub_area_normalized']),
                                           np.std(df['dis_to_hub_area_normalized']), np.quantile(df['dis_to_hub_area_normalized'].tolist(), np.arange(0, 1.1, 0.25))[2],
                                           np.mean(df['area_nuclear']), np.mean(df['mean_int_nuclear']),
                                           np.mean(df['mean_int_DNAFISH']), np.mean(df['n_ecDNA']),
                                           np.mean(df['mean_int_ecDNA']), np.mean(df['total_area_ecDNA']),
                                           np.mean(df['total_area_ratio_ecDNA'])]
    else:
        df_mean.loc[len(df_mean.index)] = [sample, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
df_mean.to_csv('%s/summary.txt' % output_dir, index=False, sep='\t')
print("DONE!")
