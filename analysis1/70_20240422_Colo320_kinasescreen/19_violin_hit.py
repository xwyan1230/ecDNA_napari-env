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

batch = 'point5uM_24hr'
n_feature_group = 2

data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
data = pd.concat([data, data_cc], axis=1)
data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]

"""features = ['n_filtered', 'n_neg', 'n_pos', 'fov_hoechst', 'pos_vs_neg',
            'per_neg_G1', 'per_neg_G1S', 'per_neg_S', 'per_neg_G2M',
            'per_pos_G1', 'per_pos_G1S', 'per_pos_S', 'per_pos_G2M']"""
features = ['n_filtered', 'n_neg', 'n_pos', 'fov_hoechst', 'pos_vs_neg']
features1 = ['per_neg_G1', 'per_neg_G1S', 'per_neg_S', 'per_neg_G2M',
            'per_pos_G1', 'per_pos_G1S', 'per_pos_S', 'per_pos_G2M']

data_flt = pd.DataFrame()
for feature in features:
    data_temp = data[(data['log2fc_%s' % feature] <= -1) | (data['log2fc_%s' % feature] >= 1)]
    data_flt = pd.concat([data_flt, data_temp], axis=0).drop_duplicates().reset_index(drop=True)

if n_feature_group == 2:
    data_flt1 = pd.DataFrame()
    for feature in features1:
        data_temp = data[(data['log2fc_%s' % feature] <= -1) | (data['log2fc_%s' % feature] >= 1)]
        data_flt1 = pd.concat([data_flt1, data_temp], axis=0).drop_duplicates().reset_index(drop=True)
        data_fltA = data_flt.merge(data_flt1, how='inner').reset_index(drop=True)
        data_fltB = data_flt[~data_flt['label'].isin(data_flt1['label'].tolist())].copy().reset_index(drop=True)

data_fltA = data_fltA.sort_values(by=['log2fc_fov_hoechst']).reset_index(drop=True)
data_fltA['identity'] = ['survival and cell cycle'] * len(data_fltA)
data_fltB = data_fltB.sort_values(by=['log2fc_fov_hoechst']).reset_index(drop=True)
data_fltB['identity'] = ['survival only'] * len(data_fltB)

df = pd.concat([data_fltA, data_fltB], axis=0).reset_index(drop=True)

features = ['n_filtered', 'n_neg', 'n_pos', 'fov_hoechst']

for feature in features:
    fig, ax = plt.subplots(figsize=(3, 5))
    fig.subplots_adjust(left=0.4)
    sns.boxplot(data=df, x='identity', y='log2fc_%s' % feature)
    sns.swarmplot(data=df, x='identity', y='log2fc_%s' % feature, color=(255/255, 140/255, 0/255))
    if not os.path.exists('%s/%s/' % (output_dir, batch)):
        os.makedirs('%s/%s/' % (output_dir, batch))
    plt.savefig('%s/%s/%s_2survivalgroup_%s.pdf' % (output_dir, batch, batch, feature))
    plt.show()

feature = 'pos_vs_neg'
fig, ax = plt.subplots(figsize=(3, 5))
fig.subplots_adjust(left=0.4)
sns.boxplot(data=df[df['mean_n_filtered'] > 200], x='identity', y='log2fc_%s' % feature)
sns.swarmplot(data=df[df['mean_n_filtered'] > 200], x='identity', y='log2fc_%s' % feature, color=(255/255, 140/255, 0/255))
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_2survivalgroup_%s.pdf' % (output_dir, batch, batch, feature))
plt.show()





