import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
# from shared.sinaplot import sinaplot
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import scipy.stats as stats
import shared.dataframe as dat
import shared.math as mat
import math
import seaborn as sns
from scipy.stats import iqr
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

folders = ['48hr_density']
# skip_lst = ['XY141', 'XY198']
skip_lst = []
density_lst = ['16k', '14k', '12k', '10k', '8k', '6k', '4k', '3k', '2k', '1k']
size_lst = [20] *10
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
fov_area = 2.336  # mm2
well_area = 11
expected_density = [16000/11, 14000/11, 12000/11, 10000/11, 8000/11, 6000/11, 4000/11, 3000/11, 2000/11, 1000/11]

data = pd.DataFrame()
for folder in folders:
    df = pd.read_csv('%s/%s/%s.txt' % (data_dir, folder, folder), na_values=['.'], sep='\t')
    df['log10_fov_hoechst'] = np.log10(df['fov_hoechst'])
    data = pd.concat([data, df], axis=0)

data['log10_fov_hoechst_per_well'] = data['log10_fov_hoechst']*well_area/fov_area
data['n_filtered_per_well'] = data['n_filtered']*well_area/fov_area

plt.subplots(figsize=(9, 7))
feature = 'log10_fov_hoechst_per_well'
for i in range(len(density_lst)):
    # plt.axhline(y=expected_density[i], color=rainboo_colors[i], linestyle='--')
    sns.scatterplot(data=data[~(data['sample'].isin(skip_lst)) & (data['density'] == density_lst[i])], x=feature, y='n_filtered_per_well', alpha=1, s=40, color= rainboo_colors[i])
# plt.xlim([0, 3])
plt.ylim([0, 18000])
plt.savefig('%s/%s_n_filtered_vs_%s.pdf' % (output_dir, folder, feature))
plt.show()






