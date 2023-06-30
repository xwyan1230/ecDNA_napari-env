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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

feature_lst = ['n', 'mean_area_nuclear', 'mean_total_area_ecDNA', 'mean_total_area_ratio_ecDNA',
               'mean_mean_int_DNAFISH', 'mean_r10', 'mean_r16', 'mean_n_ecDNA',
               'n_point2', 'mean_r10_point2', 'mean_r16_point2',
               'mean_averageD_point2', 'std_averageD_point2', 'mean_n_ecDNA_point2', 'per_n_ecDNA10_point2']
row_lst = ['C', 'D']
column_lst = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
df_mean = pd.read_csv('%ssummary_mean.txt' % data_dir, na_values=['.'], sep='\t')
df_thresh = pd.read_csv('%sthresh.txt' % data_dir, na_values=['.'], sep='\t')

for feature in feature_lst:
    feature_lst = [[] for _ in range(len(row_lst))]
    for r in range(len(row_lst)):
        for c in range(len(column_lst)):
            sample = '%s%s' % (row_lst[r], column_lst[c])
            if sample not in df_mean['sample'].tolist():
                feature_lst[r].append(-1)
            else:
                feature_lst[r].append(df_mean[df_mean['sample'] == sample][feature].tolist()[0])
    df_feature = pd.DataFrame(feature_lst)
    df_feature.index = row_lst
    df_feature.columns = ['%s' % elem for elem in column_lst]

    plt.subplots(figsize=(9, 6))
    ax1 = sns.heatmap(df_feature, cbar=0, linewidths=2, vmax=df_thresh[feature][1], vmin=df_thresh[feature][0], square=True, cmap='coolwarm', annot=True, fmt='.2f')
    if not os.path.exists("%s/heatmap/" % output_dir):
        os.makedirs("%s/heatmap/" % output_dir)
    plt.savefig('%s/heatmap/%s.pdf' % (output_dir, feature))
    plt.close()
