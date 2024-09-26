# Fig 1c

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
data_dir1 = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/processed/"
# output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/figures/"
batch = '48hr_density'

# PARAMETERS
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
fov_area = 2.336  # mm2
well_area = 11
if batch == '24hr_density':
    density_lst = [16, 14, 12, 10, 8, 6, 4, 3, 2, 1]
    ylim_val = [0, 16000]
else:
    density_lst = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    ylim_val = [0, 18000]

data = pd.read_csv('%s/01_summary/%s_normalize.txt' % (data_dir, batch), na_values=['.'], sep='\t')
data['log10_fov_hoechst_normalized_per_well'] = data['log10_fov_hoechst_normalized']*well_area/fov_area
data['n_filtered_normalized_per_well'] = data['n_filtered_normalized']*well_area/fov_area

plt.subplots(figsize=(9, 7))
feature = 'log10_fov_hoechst_normalized_per_well'
feature1 = 'n_filtered_normalized_per_well'
for i in range(len(density_lst)):
    sns.scatterplot(data=data[data['density'] == density_lst[i]], x=feature, y=feature1, alpha=1, s=40, color=rainboo_colors[i])
plt.ylim(ylim_val)
plt.show()