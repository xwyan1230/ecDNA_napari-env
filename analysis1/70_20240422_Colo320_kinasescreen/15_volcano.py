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

features = ['n_filtered', 'n_neg', 'n_pos', 'fov_hoechst', 'pos_vs_neg',
            'per_neg_G1', 'per_neg_G1S', 'per_neg_S', 'per_neg_G2M',
            'per_pos_G1', 'per_pos_G1S', 'per_pos_S', 'per_pos_G2M']
for feature in features:
    plt.subplots(figsize=(9, 9))
    sns.scatterplot(data=data, x='log2fc_%s' % feature, y='mlog10p_%s' % feature, alpha=0.5, s=30)
    sns.scatterplot(data=data_ctrl, x='log2fc_%s' % feature, y='mlog10p_%s' % feature, alpha=0.5, s=30, c='r')
    data_flt = data[(data['log2fc_%s' % feature] < -1) | (data['log2fc_%s' % feature] > 1)].copy().reset_index(drop=True)
    for i in range(len(data_flt)):
        plt.text(x=data_flt['log2fc_%s' % feature][i] + 0.03, y=data_flt['mlog10p_%s' % feature][i], s=data_flt['label'][i],
                 size=4, color=(0 / 255, 191 / 255, 255 / 255))
    if not os.path.exists('%s/%s/' % (output_dir, batch)):
        os.makedirs('%s/%s/' % (output_dir, batch))
    plt.savefig('%s/%s/%s_volcano_%s.pdf' % (output_dir, batch, batch, feature))
    plt.show()






