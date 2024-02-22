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
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

feature_lst = ['n', 'mean_area_nuclear', 'mean_total_area_ecDNA', 'mean_total_area_ratio_ecDNA',
               'mean_mean_int_DNAFISH', 'mean_r1', 'mean_r2', 'mean_r3', 'mean_e1', 'mean_e2', 'mean_e3', 'mean_n_ecDNA',
               'n_point2',
               'mean_averageD_point2', 'std_averageD_point2', 'mean_n_ecDNA_point2', 'per_n_ecDNA10_point2']

df_mean = pd.read_csv('%ssummary_mean_GFP.txt' % data_dir, na_values=['.'], sep='\t')
df_thresh = pd.read_csv('%sthresh.txt' % data_dir, na_values=['.'], sep='\t')
df_seq = pd.read_csv('%sseq_GFP.txt' % data_dir, na_values=['.'], sep='\t')

for feature in feature_lst:
    temp = []
    for i in range(len(df_seq)):
        temp.append(df_mean[df_mean['sample'] == df_seq['location'][i]][feature].tolist()[0])
    df_feature = pd.DataFrame(temp)
    df_feature.index = df_seq['gene'].tolist()
    df_feature.columns = [feature]

    fig, ax = plt.subplots(figsize=(9, 14))
    fig.subplots_adjust(left=0.3)
    ax1 = sns.heatmap(df_feature, cbar=0, linewidths=2, vmax=df_thresh[feature][1], vmin=df_thresh[feature][0], square=True, cmap='coolwarm', annot=False, fmt='.2f') # annot=True
    if not os.path.exists("%s/heatmap_line/" % output_dir):
        os.makedirs("%s/heatmap_line/" % output_dir)
    plt.savefig('%s/heatmap_line/%s_seq_GFP.pdf' % (output_dir, feature))
    plt.close()
