import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import matplotlib as mpl
import math
import numpy as np
import os
from skimage.measure import label, regionprops
from shared.sinaplot import sinaplot

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230504_analysis_DM_6hr/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

# new approach
hue_order1 = ['background', 'DNAFISH']
row_lst = ['B', 'C', 'D', 'E', 'F', 'G']
column_lst = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
sample_lst_total = []
for row in row_lst:
    sample_lst_total = sample_lst_total + ['%s%s' % (row, column) for column in column_lst]
ctrl_lst = ['C3', 'C10', 'D6', 'F2', 'F10']

# heatmap
x = np.arange(0.025, 1, 0.05)
x_label = 'relative r'
up = 20

# seg
for current_sample in sample_lst_total:
    data = pd.DataFrame()
    sample_lst = ctrl_lst + [current_sample]
    for sample in sample_lst:
        if os.path.exists('%s/txt/%s_n4.txt' % (data_dir1, sample)):
            data_temp = pd.read_csv('%s/txt/%s_n4.txt' % (data_dir1, sample), na_values=['.'], sep='\t')
            data = pd.concat([data, data_temp], axis=0)

    df = data
    df['r'] = np.sqrt(df['area_nuclear'] / math.pi)
    df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA'] / math.pi)
    df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt'] / df['r']
    df['dis_to_hub_area_normalized'] = df['dis_to_hub_area'] / df['r']

    plt.subplots(figsize=(12, len(sample_lst)))
    sinaplot(data=df, x='well', y='dis_to_hub_area_normalized', violin=False, scale='area')
    if not os.path.exists("%s/cluster/" % output_dir):
        os.makedirs("%s/cluster/" % output_dir)
    plt.savefig('%s/cluster/%s_cluster.pdf' % (output_dir, current_sample))
    plt.close()

