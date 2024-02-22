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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240207_analysis_IFandRNAFISH_test/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

folder = 'RNAFISH'
samples = ['Ctrl', 'GFP', 'mCherry', 'mix_96', 'mix_384']
# samples = ['Ctrl', 'mix_96', 'mix_384']

df = pd.DataFrame()
for sample in samples:
    df1 = pd.read_csv('%s%s/txt/%s.txt' % (data_dir, folder, sample), na_values=['.'], sep='\t')
    df1['sample'] = [sample] * len(df1)
    df = pd.concat([df, df1], axis=0)

df['log10_GFP'] = np.log10(df['GFP'])
df['log10_mCherry'] = np.log10(df['mCherry'])
df['log10_RNAFISH'] = np.log10(df['RNAFISH'])

if not os.path.exists("%s%s/" % (output_dir, folder)):
    os.makedirs("%s%s/" % (output_dir, folder))

"""plt.subplots(figsize=(12, 6))
sns.violinplot(data=df, x='sample', y='log10_RNAFISH', alpha=0.5)
plt.savefig('%s/%s/violin_log10_RNAFISH.pdf' % (output_dir, folder))
plt.show()"""

"""for sample in samples:
    plt.subplots(figsize=(9, 6))
    sns.scatterplot(data=df[df['sample'] == sample], x='log10_mCherry', y='log10_GFP', s=5, alpha=0.5, color='r')
    sns.scatterplot(data=df[df['sample'] == 'Ctrl'], x='log10_mCherry', y='log10_GFP', s=5, alpha=0.5)
    plt.xlim([2.4, 4.8])
    plt.ylim([3.2, 4.8])
    plt.savefig('%s/%s/scatter_green-red_%s.pdf' % (output_dir, folder, sample))
    plt.show()"""

x = 'mix_96'
df1 = pd.read_csv('%s%s/txt/%s.txt' % (data_dir, folder, x), na_values=['.'], sep='\t')
df1['log10_GFP'] = np.log10(df1['GFP'])
df1['log10_mCherry'] = np.log10(df1['mCherry'])
df1['log10_RNAFISH'] = np.log10(df1['RNAFISH'])
df1['log10_int'] = np.log10(df1['GFP']+df1['mCherry'])
group_lst = []
for i in range(len(df1)):
    if (df1['log10_mCherry'][i] > 2.9) & (df1['log10_GFP'][i] < 3.6):
        group_lst.append('mCherry')
    elif df1['log10_mCherry'][i] < 2.8:
        group_lst.append('GFP')
    else:
        group_lst.append('NA')
df1['group'] = group_lst

print(len(df1))
print(len(df1[df1['group'] == 'GFP']))
print(len(df1[df1['group'] == 'mCherry']))
print(len(df1[df1['group'] == 'NA']))

"""plt.subplots(figsize=(3, 6))
sns.violinplot(data=df1[df1['group'].isin(['GFP', 'mCherry'])], x='group', y='log10_RNAFISH', alpha=0.5)
plt.savefig('%s/%s/violin_log10_RNAFISH_%s.pdf' % (output_dir, folder, x))
plt.show()"""

"""plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df1, x='log10_mCherry', y='log10_GFP', s=5, alpha=0.3)
sns.scatterplot(data=df1[df1['group'] == 'GFP'], x='log10_mCherry', y='log10_GFP', s=5, alpha=0.5, color='g')
sns.scatterplot(data=df1[df1['group'] == 'mCherry'], x='log10_mCherry', y='log10_GFP', s=5, alpha=0.5, color='r')
plt.xlim([2.4, 4.8])
plt.ylim([3.2, 4.8])
# plt.savefig('%s/%s/scatter_green-red_color_%s.pdf' % (output_dir, folder, x))
plt.show()"""

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df1, x='log10_int', y='log10_RNAFISH', s=5, alpha=0.5)
sns.scatterplot(data=df1[df1['group']=='GFP'], x='log10_int', y='log10_RNAFISH', s=5, alpha=0.5, color='g')
sns.scatterplot(data=df1[df1['group']=='mCherry'], x='log10_int', y='log10_RNAFISH', s=5, alpha=0.5, color='r')
plt.xlim([3.3, 4.8])
plt.ylim([3.1, 4.6])
plt.savefig('%s/%s/int_vs_RNAFISH_%s.pdf' % (output_dir, folder, x))
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df1[df1['group']=='GFP'], x='log10_int', y='log10_RNAFISH', s=5, alpha=0.5)
plt.xlim([3.3, 4.8])
plt.ylim([3.1, 4.6])
plt.savefig('%s/%s/GFPint_vs_RNAFISH_%s.pdf' % (output_dir, folder, x))
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df1[df1['group']=='mCherry'], x='log10_int', y='log10_RNAFISH', s=5, alpha=0.5)
plt.xlim([3.3, 4.8])
plt.ylim([3.1, 4.6])
plt.savefig('%s/%s/mCherryint_vs_RNAFISH_%s.pdf' % (output_dir, folder, x))
plt.show()