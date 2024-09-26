import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
data_dir1 = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/processed/"
# output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/figures/"
batch = '5uM_24hr'

# PARAMETERS
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
fov_area = 2.336  # mm2
well_area = 11
ylim_val = [-500, 17500]
xlim_val = [40, 50]

data = pd.read_csv('%s/01_summary/%s_normalize.txt' % (data_dir, batch), na_values=['.'], sep='\t')
data['log10_fov_hoechst_per_well'] = data['log10_fov_hoechst']*well_area/fov_area
data['n_filtered_per_well'] = data['n_filtered']*well_area/fov_area
data['log10_fov_hoechst_normalized_per_well'] = data['log10_fov_hoechst_normalized']*well_area/fov_area
data['n_filtered_normalized_per_well'] = data['n_filtered_normalized']*well_area/fov_area
acquire_lst = list(set(data['acquire'].tolist()))

plt.subplots(figsize=(9, 7))
feature = 'log10_fov_hoechst_per_well'
feature1 = 'n_filtered_per_well'
for i in range(len(acquire_lst)):
    data_temp = data[data['acquire'] == acquire_lst[i]].copy().reset_index(drop=True)
    sns.scatterplot(data=data_temp, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[2*i])
plt.ylim(ylim_val)
plt.xlim(xlim_val)
plt.show()

plt.subplots(figsize=(9, 7))
feature = 'log10_fov_hoechst_normalized_per_well'
feature1 = 'n_filtered_normalized_per_well'
for i in range(len(acquire_lst)):
    data_temp = data[data['acquire'] == acquire_lst[i]].copy().reset_index(drop=True)
    sns.scatterplot(data=data_temp, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[2*i])
plt.ylim(ylim_val)
plt.xlim(xlim_val)
plt.show()
