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

feature_lst = ['n', 'averageD_mean', 'averageD_std', 'averageD_0.5',
               'area_nuclear', 'mean_int_nuclear', 'mean_int_DNAFISH', 'n_ecDNA',
               'mean_int_ecDNA', 'total_area_ecDNA', 'total_area_ratio_ecDNA']
row_lst = ['B', 'C', 'D', 'E', 'F', 'G']
column_lst = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
df = pd.read_csv('%s/summary.txt' % data_dir1, na_values=['.'], sep='\t')

for feature in feature_lst:
    feature_lst = np.array(df[feature].tolist()).reshape((len(row_lst), len(column_lst)))
    df_feature = pd.DataFrame(feature_lst)
    df_feature.index = row_lst
    df_feature.columns = ['%s' % elem for elem in column_lst]

    plt.subplots(figsize=(9, 6))
    ax1 = sns.heatmap(df_feature, cbar=0, linewidths=2, vmax=df_feature.values.max(), vmin=df_feature.values.min(), square=True, cmap='coolwarm', annot=True, fmt='.2f')
    if not os.path.exists("%s/heatmap/" % output_dir):
        os.makedirs("%s/heatmap/" % output_dir)
    plt.savefig('%s/heatmap/%s.pdf' % (output_dir, feature))
    plt.close()
