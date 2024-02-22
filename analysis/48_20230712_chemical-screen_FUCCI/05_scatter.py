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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230712_analysis_chemical-screen_FUCCI/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

feature_lst = ['G1', 'S/G2', 'G2/M', 'G1/S']
name_lst = ['G1', 'SG2', 'G2M', 'G1S']
plate_x = 'HSR_FUCCI_24hr'
plate_y = 'DM_FUCCI_24hr'

df_x = pd.read_csv('%s%s/summary_log2FC.txt' % (data_dir, plate_x), na_values=['.'], sep='\t')
df_y = pd.read_csv('%s%s/summary_log2FC.txt' % (data_dir, plate_y), na_values=['.'], sep='\t')

df_x_sort = df_x[(df_x['total'] > 200) & (df_y['total']>200)].copy().reset_index(drop=True)
df_y_sort = df_y[(df_x['total'] > 200) & (df_y['total']>200)].copy().reset_index(drop=True)

df_x_WT = df_x_sort[df_x_sort['sample'].isin(['C3', 'C10', 'D6', 'F3', 'F10'])].copy().reset_index(drop=True)
df_y_WT = df_y_sort[df_y_sort['sample'].isin(['C3', 'C10', 'D6', 'F3', 'F10'])].copy().reset_index(drop=True)

for f in range(len(feature_lst)):
    plt.subplots(figsize=(8, 8))
    plt.scatter(df_x_sort['%s_log2FC' % feature_lst[f]], df_y_sort['%s_log2FC' % feature_lst[f]], color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df_x_sort)):
        plt.text(x=df_x_sort['%s_log2FC' % feature_lst[f]][i]+0.005, y=df_y_sort['%s_log2FC' % feature_lst[f]][i]+0.005, s=df_x_sort['sample'][i], size=5)
    plt.scatter(df_x_WT['%s_log2FC' % feature_lst[f]], df_y_WT['%s_log2FC' % feature_lst[f]], color=(0.85, 0.35, 0.25), s=20)
    plt.xlabel(plate_x)
    plt.ylabel(plate_y)
    """limit1 = max(abs(df1[feature1[f]].max()), abs(df1[feature1[f]].min()))
    plt.xlim([-(limit1+0.01), (limit1+0.01)])
    limit2 = max(abs(df2[feature2[f]].max()), abs(df2[feature2[f]].min()))
    plt.ylim([-(limit2 + 0.01), (limit2 + 0.01)])"""
    limit = max(abs(df_x_sort['%s_log2FC' % feature_lst[f]].max()), abs(df_x_sort['%s_log2FC' % feature_lst[f]].min()), abs(df_y_sort['%s_log2FC' % feature_lst[f]].max()), abs(df_y_sort['%s_log2FC' % feature_lst[f]].min()))
    plt.xlim([-(limit + 0.01), (limit + 0.01)])
    plt.ylim([-(limit + 0.01), (limit + 0.01)])
    plt.savefig('%s/%s_vs_%s_%s_pw.pdf' % (output_dir, plate_x, plate_y, '%s_log2FC' % name_lst[f]))
    plt.show()