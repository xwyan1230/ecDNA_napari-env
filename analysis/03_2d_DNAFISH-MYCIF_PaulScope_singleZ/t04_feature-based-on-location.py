import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import scipy.cluster.hierarchy as shc
from pandas.plotting import parallel_coordinates
import matplotlib.pyplot as plt
import matplotlib.cm as pcm
import seaborn as sns
import shared.dataframe as dat
import shared.display as dis
import numpy as np
import os

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220726_BRDfamily_screen/"
save_folder = "%sv1_location/" % master_folder
if not os.path.exists(save_folder):
    os.makedirs(save_folder)
sample_row_lst = ['B', 'C', 'D', 'E', 'F', 'G']
sample_col_lst = np.arange(2, 12, 1)
sample_lst = []
for sample_row in sample_row_lst:
    for sample_col in sample_col_lst:
        sample_lst.append("%s%s" % (sample_row, sample_col))

data_mean = pd.read_csv(("%sv1_summary/summary_mean.txt" % master_folder), na_values=['.'], sep='\t')
data_gene = pd.read_csv("%sgene.txt" % master_folder, na_values=['.'], sep='\t')
genelist = list(set(data_gene['gene']))

feature = data_mean.columns.tolist()
feature.remove('sample')
feature.remove('gene')

for f in feature:
    value_lst = []
    value_lst_ori = data_mean[f].tolist()
    value_lst_ori.pop()
    for i in range(len(value_lst_ori)):
        if sample_lst[i] == 'D2':
            value_lst.append(0)
        value_lst.append(value_lst_ori[i])
    value_location_lst = [value_lst[i:i+10] for i in range(0, len(value_lst), 10)]
    print(value_location_lst)

    # heat map
    plt.subplots(figsize=(8, 15))
    ax1 = sns.heatmap(value_location_lst, cbar=0, linewidths=2, square=True, cmap='viridis')
    plt.savefig('%s%s.pdf' % (save_folder, f))
    plt.close()