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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240205_analysis_sp8_dWee1andiWee1/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder
line_colors = [[85/255, 103/255, 160/255], [85/255, 103/255, 160/255], [85/255, 103/255, 160/255]]
ignore = 0
# sns.set_palette(sns.color_palette(line_colors))

treatment = 'iWee1_24hr'

# samples = ['DMSO_24hr_1', 'dWee1_125nM_24hr_1', 'dWee1_250nM_24hr_1', 'dWee1_250nM_24hr_2', 'dWee1_500nM_24hr_1', 'dWee1_500nM_24hr_2', 'dWee1_1uM_24hr_1', 'dWee1_1uM_24hr_2', 'dWee1_2uM_24hr_1', 'dWee1_2uM_24hr_2']
# treatments = ['DMSO_24hr', 'dWee1_125nM_24hr', 'dWee1_250nM_24hr', 'dWee1_500nM_24hr', 'dWee1_1uM_24hr', 'dWee1_2uM_24hr']
samples = ['DMSO_24hr_1', 'iWee1_125nM_24hr_1', 'iWee1_125nM_24hr_2', 'iWee1_250nM_24hr_1', 'iWee1_250nM_24hr_2', 'iWee1_500nM_24hr_1', 'iWee1_500nM_24hr_4', 'iWee1_1uM_24hr_1', 'iWee1_1uM_24hr_2']
treatments = ['DMSO_24hr', 'iWee1_125nM_24hr', 'iWee1_250nM_24hr', 'iWee1_500nM_24hr', 'iWee1_1uM_24hr']

df = pd.DataFrame()
for sample in samples:
    df1 = pd.read_csv('%s/txt/%s_n0.txt' % (data_dir, sample), na_values=['.'], sep='\t')
    df = pd.concat([df, df1], axis=0)

df_sort = df[df['area_nuclear'] > 2000].copy().reset_index(drop=True)
df_sort1 = pd.DataFrame()
for i in range(len(treatments)):
    temp = df_sort[df_sort['sample'] == treatments[i]].copy().sort_values(by='n_ecDNA').reset_index(drop=True)
    # lowlimit = int(ignore/100 * len(temp))
    lowlimit = 0
    highlimit = int((1-ignore/100) * len(temp))
    temp = temp[lowlimit:highlimit].copy().reset_index(drop=True)
    df_sort1 = pd.concat([df_sort1, temp], axis=0)

"""_, p_val = stats.ttest_ind(df_sort[df_sort['sample'] == 'dWee1_250nM_24hr']['n_ecDNA'], df_sort[df_sort['sample'] == 'dWee1_1uM_24hr']['n_ecDNA'], equal_var=True)
print(p_val)"""

plt.subplots(figsize=(10, 9))
for i in range(len(treatments)):
    y = df_sort1[df_sort1['sample'] == treatments[i]]['n_ecDNA'].tolist()
    # Add some random "jitter" to the x-axis
    x = np.random.normal(i, 0.1, size=len(y))
    plt.plot(x, y, '.', color=(109/255, 155/255, 156/255), alpha=0.3)
sns.violinplot(data=df_sort1, x='sample', y='n_ecDNA', linecolor='black', color=(175/255, 212/255, 212/255), inner_kws=dict(box_width=12, whis_width=2), cut=0)
# plt.ylim([0, 30])
plt.savefig('%s/n_ecDNA/%s-vs-DMSO_violinplot_ignore%s.pdf' % (output_dir, treatment, ignore))
plt.show()

"""limit = 15
for sample in treatments:
    print(sample)
    print(len(df_sort[(df_sort['sample'] == sample) & (df_sort['n_ecDNA'] > limit)])*1.0/len(df_sort[df_sort['sample'] == sample]))"""



"""plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_sort, x='total_area_ecDNA_sqrt_normalized', y='dis_to_hub_area_normalized', hue='sample')
plt.show()"""

"""plt.subplots(figsize=(12, 9))
sinaplot(data=df_sort, x='sample', y='n_ecDNA', alpha=0.7, violin=False, scale='area')
if not os.path.exists("%s/n_ecDNA/" % output_dir):
    os.makedirs("%s/n_ecDNA/" % output_dir)
plt.ylim([0, 30])
plt.savefig('%s/n_ecDNA/%s-vs-DMSO_sinaplot.pdf' % (output_dir, treatment))
plt.show()"""

"""plt.subplots(figsize=(10, 9))
sns.boxplot(data=df_sort, x='sample', y='n_ecDNA', fill=False, showfliers=False)
for i in range(len(treatments)):
    y = df_sort[df_sort['sample'] == treatments[i]]['n_ecDNA'].tolist()
    # Add some random "jitter" to the x-axis
    x = np.random.normal(i, 0.1, size=len(y))
    plt.plot(x, y, 'r.', alpha=0.05)
# plt.ylim([0, 30])
plt.savefig('%s/n_ecDNA/%s-vs-DMSO_boxplot.pdf' % (output_dir, treatment))
plt.show()"""