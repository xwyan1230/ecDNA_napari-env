import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
# from shared.sinaplot import sinaplot
from scipy.stats import ks_2samp
import os
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sforgraph/" % master_folder

# samples
ctrl = pd.read_csv('%sCtrl.txt' % data_dir, na_values=['.'], sep='\t')
features = ['mean_n_ecDNA']

# ex_list = []
ex_list = ['D3', 'G10']

# load data
df_mean = pd.read_csv('%ssummary_mean.txt' % data_dir, na_values=['.'], sep='\t')
df_mean_ctrl = df_mean[df_mean['sample'].isin(ctrl['sample'].tolist())].copy().reset_index(drop=True)
df_mean_ex = df_mean[~df_mean['sample'].isin(ex_list)].copy().reset_index(drop=True)
df_mean_ex_resort = df_mean_ex.sort_values(by='mean_n_ecDNA').copy().reset_index(drop=True)
df_mean_ctrl_ex_resort = df_mean_ex_resort[(df_mean_ex_resort['sample'].isin(ctrl['sample'].tolist())) & (~df_mean_ex_resort['sample'].isin(ex_list))].copy()
df_mean_text = df_mean_ex_resort[df_mean_ex_resort['sample'] == 'G7'].copy()

if not os.path.exists("%s/n_ecDNA/" % output_dir):
    os.makedirs("%s/n_ecDNA/" % output_dir)

# scatter plot
for feature in features:
    plt.subplots(figsize=(9, 6))
    plt.scatter(list(df_mean_ex_resort.index), df_mean_ex_resort[feature], color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df_mean_text)):
        plt.text(x=list(df_mean_text.index)[i]-3.5, y=df_mean_text[feature].tolist()[i]-0.02, s=df_mean_text['gene'].tolist()[i], size=10)
    plt.scatter(list(df_mean_ctrl_ex_resort.index), df_mean_ctrl_ex_resort[feature], color=(0.85, 0.35, 0.25), s=20)
    limit = max(abs(df_mean_ex_resort[feature].max()), abs(df_mean_ex_resort[feature].min()))
    plt.ylim([-(limit+0.05), limit+0.05])
    plt.ylabel(feature)
    plt.savefig('%s/n_ecDNA/%s.pdf' % (output_dir, feature))
    plt.show()