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
from matplotlib_venn import venn2
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

folder = '5uM_24hr'
batches = [1, 2, 3]
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
colors = ['#bc4d4a', '#669daa', '#dd933c'] * 3
fov_area = 2.336  # mm2
well_area = 11
ylim_val = [-500, 17500]

data = pd.DataFrame()
for batch in batches:
    df = pd.read_csv('%s/%s/%s_%s_update.txt' % (data_dir, folder, folder, batch), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0).reset_index(drop=True)

features = ['n_G1', 'n_S', 'n_G2M', 'n_neg_G1', 'n_neg_S', 'n_neg_G2M', 'n_pos_G1', 'n_pos_S', 'n_pos_G2M']
feature_sum = ['n_filtered'] * 3 + ['n_neg'] * 3 + ['n_pos'] * 3
for i in range(len(features)):
    data['per_%s' % features[i]] = data[features[i]]/data[feature_sum[i]]


def get_average(data, data_rep1, data_rep2, feature):
    data['%s_rep1' % feature] = data_rep1[feature]
    data['%s_rep2' % feature] = data_rep2[feature]
    data['mean_%s' % feature] = (data['%s_rep1' % feature] + data['%s_rep2' % feature])/2
    return data


df = pd.DataFrame()
data_rep1 = data[data['rep']==1].copy().reset_index(drop=True).sort_values(by='group').reset_index(drop=True)
data_rep2 = data[data['rep']==2].copy().reset_index(drop=True).sort_values(by='group').reset_index(drop=True)
if data_rep1['group'].tolist() == data_rep2['group'].tolist():
    df['treatment'] = data_rep1['treatment']
    df['rep1_plate'] = data_rep1['plate']
    df['rep2_plate'] = data_rep2['plate']
    df['rep1_sample'] = data_rep1['sample']
    df['rep2_sample'] = data_rep2['sample']
    for feature in features:
        df = get_average(df, data_rep1, data_rep2, 'per_%s' % feature)
df.to_csv('%s/%s/%s_cellcycle_average.txt' % (data_dir, folder, folder), index=False, sep='\t')

data1 = data[data['n_filtered']>500].copy()
data_flt = data1[((data1['per_green']<0.1) & (data1['per_red']<0.1))|(data1['per_NA']>0.5)].copy()
print(data_flt['treatment'])
print(data_flt['rep'])
print(data_flt['plate'])
print(data_flt['sample'])

print(len(data))
data = data.drop(data_flt.index)
print(len(data))
print(len(df))
for i in data_flt.index:
    if data_flt['rep'][i]==1:
        df_temp = df[(df['rep1_plate'] == data_flt['plate'][i])&(df['rep1_sample'] == data_flt['sample'][i])]
        print(list(df_temp.index))
        df = df.drop(df_temp.index)
    else:
        df_temp = df[(df['rep2_plate'] == data_flt['plate'][i]) & (df['rep2_sample'] == data_flt['sample'][i])]
        print(list(df_temp.index))
        df = df.drop(df_temp.index)
df = df.reset_index(drop=True)
print(len(df))

for k in range(len(features)):
    f = features[k]
    rep1_std = np.std(data[(data['rep'] == 1) & (data['treatment'] == 'DMSO')]['per_%s' % f])
    rep2_std = np.std(data[(data['rep'] == 2) & (data['treatment'] == 'DMSO')]['per_%s' % f])
    rep1_mean = np.mean(data[(data['rep'] == 1) & (data['treatment'] == 'DMSO')]['per_%s' % f])
    rep2_mean = np.mean(data[(data['rep'] == 2) & (data['treatment'] == 'DMSO')]['per_%s' % f])
    print(rep1_std)
    print(rep1_mean)
    print(rep2_std)
    print(rep2_mean)

    plt.subplots(figsize=(9, 7))
    feature = 'per_%s_rep1' % f
    feature1 = 'per_%s_rep2' % f
    data_hoechst_hit = df[((df[feature] < (rep1_mean-3*rep1_std)) & (df[feature1] < (rep2_mean-3*rep2_std))) | ((df[feature] > (rep1_mean+3*rep1_std)) & (df[feature1] > (rep2_mean+3*rep2_std)))].copy().reset_index(drop=True)
    data_hoechst_hit.to_csv('%s/%s/%s_%s_hit.txt' % (data_dir, folder, folder, feature), index=False, sep='\t')
    sns.scatterplot(data=df, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
    sns.scatterplot(data=df[df['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
    sns.scatterplot(data=data_hoechst_hit, x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[1])
    plt.axvline(x=rep1_mean-3*rep1_std, linestyle='--', color=rainboo_colors[0])
    plt.axhline(y=rep2_mean-3*rep2_std, linestyle='--', color=rainboo_colors[0])
    plt.axvline(x=rep1_mean + 3 * rep1_std, linestyle='--', color=rainboo_colors[0])
    plt.axhline(y=rep2_mean + 3 * rep2_std, linestyle='--', color=rainboo_colors[0])
    data_flt = data_hoechst_hit
    for i in range(len(data_flt)):
        plt.text(x=data_flt[feature][i] + 0.01, y=data_flt[feature1][i], s=data_flt['treatment'][i],
                 size=6, color=rainboo_colors[1])
    # plt.ylim([40, 50])
    # plt.xlim([40, 50])
    if not os.path.exists('%s/%s/cellcycle/' % (output_dir, folder)):
        os.makedirs('%s/%s/cellcycle/' % (output_dir, folder))
    plt.savefig('%s/%s/cellcycle/%s_%s_rep.pdf' % (output_dir, folder, folder, feature))
    plt.show()








