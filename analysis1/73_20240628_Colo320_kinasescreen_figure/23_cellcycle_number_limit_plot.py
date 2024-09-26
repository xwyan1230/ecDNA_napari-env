import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import scipy.cluster.hierarchy as shc
from pandas.plotting import parallel_coordinates
import matplotlib.pyplot as plt
import matplotlib.cm as pcm
import matplotlib
import seaborn as sns
import shared.dataframe as dat
import shared.display as dis
import numpy as np
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

df = pd.read_csv('%s/cell-cycle_sampling.txt' % (data_dir), na_values=['.'], sep='\t')

plt.subplots(figsize=(9, 7))
# plt.bar(df[df['n_filtered'] > 500]['group'], df[df['n_filtered'] > 500]['percentage'], 0.5)
sns.violinplot(data=df, x='sample_n', y='per_G2M', color='#e3e3e3', edgecolor='#4d4e4e')
# sns.swarmplot(data=df, x='sample_n', y='per_G1', color='#bc4d4a')
plt.ylim([0, 1.05])
"""if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/cutoff_qc.pdf' % output_dir)"""
plt.show()