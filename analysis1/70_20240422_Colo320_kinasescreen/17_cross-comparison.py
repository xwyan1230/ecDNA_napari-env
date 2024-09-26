import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
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
    data_flt = data_flt.merge(data_flt1, how='inner').reset_index(drop=True)

data_flt = data_flt.sort_values(by=['log2fc_fov_hoechst']).reset_index(drop=True)

data_hit_48 = data_flt

batch1 = 'point5uM_24hr'

data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch1, batch1), na_values=['.'], sep='\t')
data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch1, batch1), na_values=['.'], sep='\t')
data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
data = pd.concat([data, data_cc], axis=1)
data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]


data_flt = pd.DataFrame()
for feature in features:
    data_temp = data[(data['log2fc_%s' % feature] <= -1) | (data['log2fc_%s' % feature] >= 1)]
    data_flt = pd.concat([data_flt, data_temp], axis=0).drop_duplicates().reset_index(drop=True)

if n_feature_group == 2:
    data_flt1 = pd.DataFrame()
    for feature in features1:
        data_temp = data[(data['log2fc_%s' % feature] <= -1) | (data['log2fc_%s' % feature] >= 1)]
        data_flt1 = pd.concat([data_flt1, data_temp], axis=0).drop_duplicates().reset_index(drop=True)
    data_flt = data_flt.merge(data_flt1, how='inner').reset_index(drop=True)

data_flt = data_flt.sort_values(by=['log2fc_fov_hoechst']).reset_index(drop=True)

data_hit_24 = data_flt

setA = data_hit_48['label'].tolist()
setB = data_hit_24['label'].tolist()

red = []  # in A not in B
orange = []  # common
green = []  # in B not in A

for i in setA:
    if i in setB:
        orange.append(i)
    else:
        red.append(i)

for j in setB:
    if j not in setA:
        green.append(j)

print('red')
print(red)
print('orange')
print(orange)
print('green')
print(green)

plt.subplots(figsize=(9, 9))
venn2([set(setA), set(setB)])
plt.savefig('%s/%s_vs_%s_survival_and_cellcycle.pdf' % (output_dir, batch, batch1))
plt.show()
