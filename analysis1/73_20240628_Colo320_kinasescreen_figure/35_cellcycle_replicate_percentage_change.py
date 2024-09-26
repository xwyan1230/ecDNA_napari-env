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
exclude_index = [154, 188]
exclude_group = [155, 189]
fov_area = 2.336  # mm2
well_area = 11
ylim_val = [-500, 17500]

data = pd.DataFrame()
for batch in batches:
    df = pd.read_csv('%s/%s/%s_%s_update1.txt' % (data_dir, folder, folder, batch), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0).reset_index(drop=True)

features = ['n_G1', 'n_G1S', 'n_S', 'n_G2M', 'n_G2MG1', 'n_neg_G1', 'n_neg_G1S', 'n_neg_S', 'n_neg_G2M', 'n_neg_G2MG1',
            'n_pos_G1', 'n_pos_G1S', 'n_pos_S', 'n_pos_G2M', 'n_pos_G2MG1']
feature_sum = ['n_filtered'] * 5 + ['n_neg'] * 5 + ['n_pos'] * 5

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
    df['index'] = df.index
    df['rep1_plate'] = data_rep1['plate']
    df['rep2_plate'] = data_rep2['plate']
    df['rep1_sample'] = data_rep1['sample']
    df['rep2_sample'] = data_rep2['sample']
    for feature in features:
        df = get_average(df, data_rep1, data_rep2, 'per_%s' % feature)

print(len(data))
data_drop = data.copy()
data_drop = data_drop[~data_drop['group'].isin(exclude_group)].copy().reset_index(drop=True)
print(len(data_drop))

for k in range(len(features)):
    f = features[k]
    rep1_std = np.std(data_drop[(data_drop['rep'] == 1) & (data_drop['treatment'] == 'DMSO')]['per_%s' % f])
    rep2_std = np.std(data_drop[(data_drop['rep'] == 2) & (data_drop['treatment'] == 'DMSO')]['per_%s' % f])
    rep1_mean = np.mean(data_drop[(data_drop['rep'] == 1) & (data_drop['treatment'] == 'DMSO')]['per_%s' % f])
    rep2_mean = np.mean(data_drop[(data_drop['rep'] == 2) & (data_drop['treatment'] == 'DMSO')]['per_%s' % f])
    print(rep1_std)
    print(rep1_mean)
    print(rep2_std)
    print(rep2_mean)
    df['delta_%s_rep1' % f] = [df['per_%s_rep1' % f][i] - rep1_mean for i in range(len(df))]
    df['delta_%s_rep2' % f] = [df['per_%s_rep2' % f][i] - rep2_mean for i in range(len(df))]
    df['mean_delta_%s' % f] = (df['delta_%s_rep1' % f] + df['delta_%s_rep2' % f])/2
    """df['perchange_%s_rep1' % f] = [(df['per_%s_rep1' % f][i] - rep1_mean)/rep1_mean if df['per_%s_rep1' % f][i] <= rep1_mean
                                   else (df['per_%s_rep1' % f][i] - rep1_mean)/(1-rep1_mean) for i in range(len(df))]
    df['perchange_%s_rep2' % f] = [
        (df['per_%s_rep2' % f][i] - rep2_mean) / rep2_mean if df['per_%s_rep2' % f][i] <= rep2_mean
        else (df['per_%s_rep2' % f][i] - rep2_mean) / (1 - rep2_mean) for i in range(len(df))]
    df['mean_perchange_%s' % f] = (df['perchange_%s_rep1' % f] + df['perchange_%s_rep2' % f])/2"""

    print(len(df))
    df_drop = df.copy()
    df_drop = df_drop.drop(exclude_index)
    print(len(df_drop))

    plt.subplots(figsize=(9, 7))
    feature = 'per_%s_rep1' % f
    feature1 = 'per_%s_rep2' % f
    # data_hoechst_hit = df_drop[((df_drop[feature] < (rep1_mean-3*rep1_std)) & (df_drop[feature1] < (rep2_mean-3*rep2_std))) | ((df_drop[feature] > (rep1_mean+3*rep1_std)) & (df_drop[feature1] > (rep2_mean+3*rep2_std)))].copy().reset_index(drop=True)
    # data_hoechst_hit.to_csv('%s/%s/%s_%s_hit.txt' % (data_dir, folder, folder, feature), index=False, sep='\t')
    sns.scatterplot(data=df_drop, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
    sns.scatterplot(data=df_drop[df_drop['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
    # sns.scatterplot(data=data_hoechst_hit, x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[1])
    # plt.axvline(x=rep1_mean-3*rep1_std, linestyle='--', color=rainboo_colors[0])
    # plt.axhline(y=rep2_mean-3*rep2_std, linestyle='--', color=rainboo_colors[0])
    # plt.axvline(x=rep1_mean + 3 * rep1_std, linestyle='--', color=rainboo_colors[0])
    # plt.axhline(y=rep2_mean + 3 * rep2_std, linestyle='--', color=rainboo_colors[0])
    data_flt = df_drop[df_drop[feature]>0.2].copy().reset_index(drop=True)
    for i in range(len(data_flt)):
        plt.text(x=data_flt[feature][i] + 0.01, y=data_flt[feature1][i], s=data_flt['treatment'][i],
                 size=6, color=rainboo_colors[1])
    plt.ylim([0, 1])
    plt.xlim([0, 1])
    if not os.path.exists('%s/%s/cellcycle/' % (output_dir, folder)):
        os.makedirs('%s/%s/cellcycle/' % (output_dir, folder))
    plt.savefig('%s/%s/cellcycle/%s_%s_rep_update1.pdf' % (output_dir, folder, folder, feature))
    plt.show()

    """plt.subplots(figsize=(9, 7))
    feature = 'per_%s_rep1' % f
    feature1 = 'per_%s_rep2' % f
    sns.scatterplot(data=df_drop, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
    sns.scatterplot(data=df_drop[df_drop['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20,
                    color=rainboo_colors[6])
    plt.ylim([0, 1])
    plt.xlim([0, 1])
    if not os.path.exists('%s/%s/cellcycle/' % (output_dir, folder)):
        os.makedirs('%s/%s/cellcycle/' % (output_dir, folder))
    plt.savefig('%s/%s/cellcycle/%s_%s_rep_woheat.pdf' % (output_dir, folder, folder, feature))
    plt.show()"""

    """plt.subplots(figsize=(9, 7))
    feature = 'perchange_%s_rep1' % f
    feature1 = 'perchange_%s_rep2' % f
    sns.scatterplot(data=df_drop, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
    sns.scatterplot(data=df_drop[df_drop['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
    # plt.ylim([40, 50])
    # plt.xlim([40, 50])
    if not os.path.exists('%s/%s/cellcycle/' % (output_dir, folder)):
        os.makedirs('%s/%s/cellcycle/' % (output_dir, folder))
    plt.savefig('%s/%s/cellcycle/%s_%s_rep_perchange.pdf' % (output_dir, folder, folder, feature))
    plt.show()"""

df.to_csv('%s/%s/%s_cellcycle_average_update1.txt' % (data_dir, folder, folder), index=False, sep='\t')








