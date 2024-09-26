import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import seaborn as sns
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

batch = '5uM_24hr'

data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
data = pd.concat([data, data_cc], axis=1)
data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]
data_ctrl = data[data['treatment'] == 'DMSO'].copy().reset_index(drop=True)
data_flt1 = data[(data['rep_neg'] == 2) & (data['rep_pos'] == 2)
                & (data['mean_n_neg_larger_than_100'] == 1) & (data['mean_n_neg_larger_than_100'] == 1)].copy().reset_index(drop=True)

features = ['G1', 'G1S', 'S', 'G2M']
for feature in features:
    plt.subplots(figsize=(9, 9))
    sns.scatterplot(data=data, x='log2fc_per_neg_%s' % feature, y='log2fc_per_pos_%s' % feature, alpha=0.3, s=30, c=(0.1, 0.1, 0.1))
    sns.scatterplot(data=data_flt1, x='log2fc_per_neg_%s' % feature, y='log2fc_per_pos_%s' % feature, alpha=0.5, s=30)
    sns.scatterplot(data=data_ctrl, x='log2fc_per_neg_%s' % feature, y='log2fc_per_pos_%s' % feature, s=30, c='r')
    plt.xlim([-16, 3])
    plt.ylim([-16, 3])
    data_flt = data_flt1[(data_flt1['log2fc_per_neg_%s' % feature] <= -1) | (data_flt1['log2fc_per_neg_%s' % feature] >= 1)
                    | (data_flt1['log2fc_per_pos_%s' % feature] <= -1) | (data_flt1['log2fc_per_pos_%s' % feature] >= 1)].copy().reset_index(drop=True)
    for i in range(len(data_flt)):
        plt.text(x=data_flt['log2fc_per_neg_%s' % feature][i] + 0.03, y=data_flt['log2fc_per_pos_%s' % feature][i], s=data_flt['label'][i],
                 size=4, color=(0 / 255, 191 / 255, 255 / 255))
    if not os.path.exists('%s/%s/' % (output_dir, batch)):
        os.makedirs('%s/%s/' % (output_dir, batch))
    plt.savefig('%s/%s/%s_cellcycle_%s.pdf' % (output_dir, batch, batch, feature))
    plt.show()




