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

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220526_flowFISH_topHits_screen/"
save_folder = "%sv5_location/" % master_folder
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

data_mean = pd.read_csv(("%sv5_summary/summary_mean.txt" % master_folder), na_values=['.'], sep='\t')
data_gene = pd.read_csv("%sgene.txt" % master_folder, na_values=['.'], sep='\t')
genelist = list(set(data_gene['gene']))

feature = data_mean.columns.tolist()
feature.remove('sample')
feature.remove('gene')

for f in feature:
    value_lst = data_mean[f].tolist()
    value_lst.pop()
    value_location_lst = [value_lst[i:i+8] for i in range(0, len(value_lst), 8)]
    print(value_location_lst)

    # heat map
    plt.subplots(figsize=(8, 15))
    ax1 = sns.heatmap(value_location_lst, cbar=0, linewidths=2, square=True, cmap='viridis')
    plt.savefig('%s%s.pdf' % (save_folder, f))
    plt.close()