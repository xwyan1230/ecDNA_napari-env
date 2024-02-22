import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
from shared.sinaplot import sinaplot
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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240205_analysis_sp8_dWee1andiWee1/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

samples = ['DMSO_24hr_1', 'dWee1_1uM_24hr_1', 'dWee1_1uM_24hr_2']

df = pd.DataFrame()
for sample in samples:
    df1 = pd.read_csv('%s/txt/%s_n0.txt' % (data_dir, sample), na_values=['.'], sep='\t')
    df = pd.concat([df, df1], axis=0)

df['r'] = np.sqrt(df['area_nuclear'] / math.pi)
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA'] / math.pi)
df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt'] / df['r']
df['dis_to_hub_area_normalized'] = df['dis_to_hub_area'] / df['r']

df_sort = df.copy()

plt.subplots(figsize=(4, 9))
sinaplot(data=df_sort, x='sample', y='n_ecDNA', alpha=0.7, violin=False, scale='area')
plt.show()

plt.subplots(figsize=(4, 9))
sns.violinplot(data=df_sort, x='sample', y='n_ecDNA', cut=0)
plt.show()

print(len(df_sort[df_sort['sample'] == 'DMSO_24hr']))
print(len(df_sort[df_sort['sample'] == 'dWee1_1uM_24hr']))
print(stats.ttest_ind(df_sort[df_sort['sample'] == 'dWee1_1uM_24hr']['n_ecDNA'], df_sort[df_sort['sample'] == 'DMSO_24hr']['n_ecDNA'], equal_var=True))
print(stats.mannwhitneyu(df_sort[df_sort['sample'] == 'dWee1_1uM_24hr']['n_ecDNA'], df_sort[df_sort['sample'] == 'DMSO_24hr']['n_ecDNA'], use_continuity=False))
# if not os.path.exists("%s/cluster/" % output_dir):
#     os.makedirs("%s/cluster/" % output_dir)
# plt.savefig('%s/cluster/%s_cluster.pdf' % (output_dir, current_sample))

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_sort, x='total_area_ecDNA_sqrt_normalized', y='dis_to_hub_area_normalized', hue='sample')
plt.show()