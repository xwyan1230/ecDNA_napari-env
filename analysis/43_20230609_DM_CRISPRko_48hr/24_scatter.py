import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
from shared.sinaplot import sinaplot
from scipy.stats import ks_2samp
import os
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

row_lst = ['C', 'D', 'E', 'F', 'G']
column_lst = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

plt.subplots(figsize=(9, 6))

for r in range(len(row_lst)):
    for c in range(len(column_lst)):
        sample = '%s%s' % (row_lst[r], column_lst[c])
        if os.path.exists('%s/%s/%s_summary.txt' % (data_dir, sample, sample)):
            df = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
            plt.scatter(df['mean_mean_int_DNAFISH'][0], df['mean_averageD_point2'][0], color=(220/255, 20/255, 60/255), s=20)
            plt.text(x=df['mean_mean_int_DNAFISH'][0] + 5, y=df['mean_averageD_point2'][0] + 0.01, s='%s_mCherry' % sample, size=5)
            plt.scatter(df['mean_mean_int_DNAFISH'][1], df['mean_averageD_point2'][1], color=(154 / 255, 205 / 255, 50 / 255), s=20)
            plt.text(x=df['mean_mean_int_DNAFISH'][1] + 5, y=df['mean_averageD_point2'][1] + 0.01, s='%s_GFP' % sample, size=5)
plt.xlabel('mean_int_DNAFISH')
plt.ylabel('mean_averageD_point2')
plt.savefig('%s/volcano/scatter.pdf' % output_dir)
plt.show()