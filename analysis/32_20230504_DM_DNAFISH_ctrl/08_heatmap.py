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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230504_analysis_DM_DNAFISH_ctrl/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

exp = '20230409_DM_plate_DNAFISH_control'

feature_lst = ['n', 'averageD_mean', 'averageD_std', 'averageD_0.1', 'averageD_0.25', 'averageD_0.5',
               'area_nuclear', 'mean_int_nuclear', 'mean_int_DNAFISH', 'n_ecDNA',
               'mean_int_ecDNA', 'total_area_ecDNA', 'total_area_ratio_ecDNA']
row_lst = ['B', 'C', 'D', 'E', 'F', 'G']
column_lst = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
"""sample_lst = []
for row in row_lst:
    sample_lst = sample_lst + ['%s%s' % (row, column) for column in column_lst]"""

for feature in feature_lst:
    feature_lst = [[] for _ in range(len(row_lst))]
    for r in range(len(row_lst)):
        for c in range(len(column_lst)):
            sample = '%s%s' % (row_lst[r], column_lst[c])
            if not os.path.exists('%s/txt/%s_n4.txt' % (data_dir1, sample)):
                feature_lst[r].append(-1)
            else:
                df = pd.read_csv('%s/txt/%s_n4.txt' % (data_dir1, sample), na_values=['.'], sep='\t')
                if feature == 'n':
                    feature_lst[r].append(len(df))
                elif 'averageD' in feature:
                    df['r'] = np.sqrt(df['area_nuclear'] / math.pi)
                    df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA'] / math.pi)
                    df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt'] / df['r']
                    df['dis_to_hub_area_normalized'] = df['dis_to_hub_area'] / df['r']
                    if feature == 'averageD_mean':
                        feature_lst[r].append(np.mean(df['dis_to_hub_area_normalized']))
                    elif feature == 'averageD_std':
                        feature_lst[r].append(np.std(df['dis_to_hub_area_normalized']))
                    elif feature == 'averageD_0.1':
                        feature_lst[r].append(np.quantile(df['dis_to_hub_area_normalized'].tolist(), np.arange(0, 1.1, 0.1))[1])
                    elif feature == 'averageD_0.25':
                        feature_lst[r].append(np.quantile(df['dis_to_hub_area_normalized'].tolist(), np.arange(0, 1.1, 0.25))[1])
                    elif feature == 'averageD_0.5':
                        feature_lst[r].append(np.quantile(df['dis_to_hub_area_normalized'].tolist(), np.arange(0, 1.1, 0.25))[2])
                else:
                    feature_lst[r].append(np.mean(df[feature]))
    df_feature = pd.DataFrame(feature_lst)
    df_feature.index = row_lst
    df_feature.columns = ['%s' % elem for elem in column_lst]

    plt.subplots(figsize=(9, 6))
    ax1 = sns.heatmap(df_feature, cbar=0, linewidths=2, vmax=df_feature.values.max(), vmin=df_feature.values.min(), square=True, cmap='coolwarm', annot=True, fmt='.2f')
    if not os.path.exists("%s/heatmap/" % output_dir):
        os.makedirs("%s/heatmap/" % output_dir)
    plt.savefig('%s/heatmap/%s.pdf' % (output_dir, feature))
    plt.close()
