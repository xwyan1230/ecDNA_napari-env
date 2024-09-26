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
import os
from matplotlib_venn import venn2

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

folder = 'point5uM_48hr'
batches = [1, 2, 3]
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
fov_area = 2.336  # mm2
well_area = 11
test3 = ['XY0%s' % (x+1) for x in range(9)] + ['XY%s' % (x+10) for x in range(5)]
ylim_val = [-500, 17500]

data = pd.DataFrame()
for batch in batches:
    df = pd.read_csv('%s/%s/%s_%s_update1.txt' % (data_dir, folder, folder, batch), na_values=['.'], sep='\t')
    df['log10_fov_hoechst'] = np.log10(df['fov_hoechst'])
    data = pd.concat([data, df], axis=0)

data['log10_fov_hoechst_per_well'] = data['log10_fov_hoechst']*well_area/fov_area
data['n_filtered_per_well'] = data['n_filtered']*well_area/fov_area

data1 = data[(data['rep'] != 0) & (data['plate'] == batches[0])].copy().reset_index(drop=True)
data2 = data[(data['rep'] != 0) & (data['plate'] == batches[1])].copy().reset_index(drop=True)
data3 = data[(data['rep'] != 0) & (data['plate'] == batches[2]) & (data['sample'].isin(test3))].copy().reset_index(drop=True)
data4 = data[(data['rep'] != 0) & (data['plate'] == batches[2]) & (~data['sample'].isin(test3))].copy().reset_index(drop=True)
data1_ctrl = data1[data1['treatment']=='DMSO'].copy().reset_index(drop=True)
data2_ctrl = data2[data2['treatment']=='DMSO'].copy().reset_index(drop=True)
data3_ctrl = data3[data3['treatment']=='DMSO'].copy().reset_index(drop=True)
data4_ctrl = data4[data4['treatment']=='DMSO'].copy().reset_index(drop=True)

"""plt.subplots(figsize=(9, 7))
feature = 'log10_fov_hoechst_per_well'
feature1 = 'n_filtered_per_well'
sns.scatterplot(data=data1, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[0])
sns.scatterplot(data=data2, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[2])
sns.scatterplot(data=data3, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[4])
sns.scatterplot(data=data4, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[6])
plt.ylim(ylim_val)
plt.xlim([40, 50])
if not os.path.exists('%s/%s/normalization/' % (output_dir, folder)):
    os.makedirs('%s/%s/normalization/' % (output_dir, folder))
plt.savefig('%s/%s/normalization/%s_before.pdf' % (output_dir, folder, folder))
plt.show()"""


def get_normalization(data, ctrl, ref):
    mean_hoechst = np.mean(ref[ref['n_filtered']>1000]['log10_fov_hoechst'])
    mean_n = np.mean(ref[ref['n_filtered']>1000]['n_filtered'])
    hoechst_ratio = np.mean(ctrl[ctrl['n_filtered']>1000]['log10_fov_hoechst'])/mean_hoechst
    n_ratio = np.mean(ctrl[ctrl['n_filtered']>1000]['n_filtered'])/mean_n
    data['log10_fov_hoechst_normalized'] = data['log10_fov_hoechst']/hoechst_ratio
    data['n_filtered_normalized'] = data['n_filtered']/n_ratio
    data['n_pos_normalized'] = data['n_pos'] / n_ratio
    data['n_neg_normalized'] = data['n_neg'] / n_ratio
    data['log10_fov_hoechst_normalized_per_well'] = data['log10_fov_hoechst_normalized'] * well_area / fov_area
    data['n_filtered_normalized_per_well'] = data['n_filtered_normalized'] * well_area / fov_area
    data['n_pos_normalized_per_well'] = data['n_pos_normalized'] * well_area / fov_area
    data['n_neg_normalized_per_well'] = data['n_neg_normalized'] * well_area / fov_area
    return data


data1 = get_normalization(data1, data1_ctrl, data1_ctrl)
data2 = get_normalization(data2, data2_ctrl, data1_ctrl)
data3 = get_normalization(data3, data3_ctrl, data1_ctrl)
data4 = get_normalization(data4, data4_ctrl, data1_ctrl)
data1_ctrl = data1[data1['treatment']=='DMSO'].copy().reset_index(drop=True)
data2_ctrl = data2[data2['treatment']=='DMSO'].copy().reset_index(drop=True)
data3_ctrl = data3[data3['treatment']=='DMSO'].copy().reset_index(drop=True)
data4_ctrl = data4[data4['treatment']=='DMSO'].copy().reset_index(drop=True)

data_all = pd.concat([data1, data2, data3, data4], axis=0)
data_all.to_csv('%s/%s/%s_normalized_update1.txt' % (data_dir, folder, folder), index=False, sep='\t')

"""plt.subplots(figsize=(9, 7))
feature = 'log10_fov_hoechst_normalized_per_well'
feature1 = 'n_filtered_normalized_per_well'
sns.scatterplot(data=data1, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[0])
sns.scatterplot(data=data2, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[2])
sns.scatterplot(data=data3, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[4])
sns.scatterplot(data=data4, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[6])
plt.ylim(ylim_val)
plt.xlim([40, 50])
plt.savefig('%s/%s/normalization/%s_after.pdf' % (output_dir, folder, folder))
plt.show()"""

"""plt.subplots(figsize=(9, 7))
feature = 'log10_fov_hoechst_normalized_per_well'
feature1 = 'n_filtered_normalized_per_well'
sns.scatterplot(data=data_all, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
sns.scatterplot(data=data_all[data_all['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
plt.ylim([-500, 16000])
plt.xlim([40, 50])
plt.savefig('%s/%s/normalization/%s_after_DMSO.pdf' % (output_dir, folder, folder))
plt.show()"""

df = pd.DataFrame()
data_rep1 = data_all[data_all['rep']==1].copy().reset_index(drop=True).sort_values(by='group').reset_index(drop=True)
data_rep2 = data_all[data_all['rep']==2].copy().reset_index(drop=True).sort_values(by='group').reset_index(drop=True)
if data_rep1['group'].tolist() == data_rep2['group'].tolist():
    df['treatment'] = data_rep1['treatment']
    df['hoechst_rep1'] = data_rep1['log10_fov_hoechst_normalized_per_well']
    df['hoechst_rep2'] = data_rep2['log10_fov_hoechst_normalized_per_well']
    df['n_rep1'] = data_rep1['n_filtered_normalized_per_well']
    df['n_rep2'] = data_rep2['n_filtered_normalized_per_well']
    df['n_pos_rep1'] = data_rep1['n_pos_normalized_per_well']
    df['n_pos_rep2'] = data_rep2['n_pos_normalized_per_well']
    df['n_neg_rep1'] = data_rep1['n_neg_normalized_per_well']
    df['n_neg_rep2'] = data_rep2['n_neg_normalized_per_well']
    df['mean_hoechst'] = (df['hoechst_rep1']+df['hoechst_rep2'])/2
    df['mean_n'] = (df['n_rep1']+df['n_rep2'])/2
    df['mean_n_pos'] = (df['n_pos_rep1']+df['n_pos_rep2'])/2
    df['mean_n_neg'] = (df['n_neg_rep1'] + df['n_neg_rep2']) / 2

features = ['log10_fov_hoechst_normalized_per_well', 'n_filtered_normalized_per_well', 'n_pos_normalized_per_well', 'n_neg_normalized_per_well']
feature_names = ['hoechst', 'n', 'n_pos', 'n_neg']
for i in range(len(features)):
    feature = features[i]
    print(feature)
    feature_name = feature_names[i]
    rep1_std = np.std(data_all[(data_all['rep'] == 1) & (data_all['treatment'] == 'DMSO')][feature])
    rep2_std = np.std(data_all[(data_all['rep'] == 2) & (data_all['treatment'] == 'DMSO')][feature])
    rep1_mean = np.mean(data_all[(data_all['rep'] == 1) & (data_all['treatment'] == 'DMSO')][feature])
    rep2_mean = np.mean(data_all[(data_all['rep'] == 2) & (data_all['treatment'] == 'DMSO')][feature])
    df['log2_%s_rep1' % feature_name] = np.log2(df['%s_rep1' % feature_name] / rep1_mean)
    df['log2_%s_rep2' % feature_name] = np.log2(df['%s_rep2' % feature_name] / rep2_mean)
    df['log2_mean_%s' % feature_name] = (df['log2_%s_rep1' % feature_name] + df['log2_%s_rep2' % feature_name]) / 2
    mean_val = np.mean(df[df['treatment'] == 'DMSO']['log2_mean_%s' % feature_name])
    std_val = np.std(df[df['treatment'] == 'DMSO']['log2_mean_%s' % feature_name])
    print(mean_val)
    print(std_val)

    plt.subplots(figsize=(9, 7))
    feature = '%s_rep1' % feature_name
    feature1 = '%s_rep2' % feature_name
    data_n_hit = df[
        (df[feature] < (rep1_mean - 3 * rep1_std)) & (df[feature1] < (rep2_mean - 3 * rep2_std))].copy().reset_index(
        drop=True)
    # data_n_hit.to_csv('%s/%s/%s_n_hit.txt' % (data_dir, folder, folder), index=False, sep='\t')
    sns.scatterplot(data=df, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
    sns.scatterplot(data=df[df['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
    sns.scatterplot(data=data_n_hit, x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[2])
    plt.axvline(x=rep1_mean - 3 * rep1_std, linestyle='--', color=rainboo_colors[0])
    plt.axhline(y=rep2_mean - 3 * rep2_std, linestyle='--', color=rainboo_colors[0])
    # plt.ylim([-500, 16000])
    # plt.xlim([-500, 16000])
    plt.savefig('%s/%s/normalization/%s_%s_rep_update1.pdf' % (output_dir, folder, folder, feature_name))
    plt.show()

    plt.subplots(figsize=(9, 7))
    feature = 'log2_%s_rep1' % feature_name
    feature1 = 'log2_%s_rep2' % feature_name
    sns.scatterplot(data=df, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
    sns.scatterplot(data=df[df['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
    # plt.ylim([40, 50])
    # plt.xlim([40, 50])
    plt.savefig('%s/%s/normalization/%s_%s_rep_log2_update1.pdf' % (output_dir, folder, folder, feature_name))
    plt.show()

df['ratio'] = df['mean_n_pos']/df['mean_n_neg']
ratio_mean = np.mean(df[df['treatment'] == 'DMSO']['ratio'])
df['log2_ratio'] = np.log2(df['ratio']/ratio_mean)
mean_val = np.mean(df[df['treatment'] == 'DMSO']['log2_ratio'])
std_val = np.std(df[df['treatment'] == 'DMSO']['log2_ratio'])
print(mean_val)
print(std_val)

df.to_csv('%s/%s/%s_average_update1.txt' % (data_dir, folder, folder), index=False, sep='\t')

plt.subplots(figsize=(9, 7))
feature = 'mean_hoechst'
feature1 = 'mean_n'
sns.scatterplot(data=df, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
sns.scatterplot(data=df[df['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
data_flt = df[df[feature]/df[feature1]<1.3].copy().reset_index(drop=True)
sns.scatterplot(data=data_flt, x=feature, y=feature1, alpha=1, s=20, color='#edb08e')  # '#d84480' '#8f44d8'
for i in range(len(data_flt)):
    plt.text(x=data_flt[feature][i] + 0.1, y=data_flt[feature1][i], s=data_flt['treatment'][i],
             size=6, color=rainboo_colors[2])
# plt.ylim([-6, 1])
# plt.xlim([-6, 1])
# plt.savefig('%s/%s/normalization/%s_n-pos_vs_n-neg.pdf' % (output_dir, folder, folder))
plt.show()

plt.subplots(figsize=(9, 7))
feature = 'log2_mean_hoechst'
feature1 = 'log2_mean_n'
sns.scatterplot(data=df, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
sns.scatterplot(data=df[df['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
data_flt = df[(df[feature1]/df[feature]>40) & (df[feature1]<-0.5)].copy().reset_index(drop=True)
sns.scatterplot(data=data_flt, x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[2])  # '#edb08e' '#d84480' '#8f44d8'
for i in range(len(data_flt)):
    plt.text(x=data_flt[feature][i] + 0.0003, y=data_flt[feature1][i], s=data_flt['treatment'][i],
             size=6, color=rainboo_colors[2])
# plt.ylim([-6, 1])
# plt.xlim([-6, 1])
plt.savefig('%s/%s/normalization/%s_hoechst_vs_n_log2.pdf' % (output_dir, folder, folder))
plt.show()

plt.subplots(figsize=(9, 7))
feature = 'mean_n_neg'
feature1 = 'mean_n_pos'
sns.scatterplot(data=df, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
sns.scatterplot(data=df[df['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
data_flt = df[df[feature]/df[feature1]<1.3].copy().reset_index(drop=True)
sns.scatterplot(data=data_flt, x=feature, y=feature1, alpha=1, s=20, color='#edb08e')  # '#d84480' '#8f44d8'
for i in range(len(data_flt)):
    plt.text(x=data_flt[feature][i] + 0.1, y=data_flt[feature1][i], s=data_flt['treatment'][i],
             size=6, color=rainboo_colors[2])
# plt.ylim([-6, 1])
# plt.xlim([-6, 1])
plt.savefig('%s/%s/normalization/%s_n-pos_vs_n-neg.pdf' % (output_dir, folder, folder))
plt.show()

plt.subplots(figsize=(9, 7))
feature = 'ratio'
feature1 = 'mean_n'
sns.scatterplot(data=df, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
sns.scatterplot(data=df[df['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
data_flt = df[df[feature]>0].copy().reset_index(drop=True)
sns.scatterplot(data=data_flt, x=feature, y=feature1, alpha=1, s=20, color='#edb08e')  # '#d84480' '#8f44d8'
for i in range(len(data_flt)):
    plt.text(x=data_flt[feature][i] + 0.1, y=data_flt[feature1][i], s=data_flt['treatment'][i],
             size=6, color=rainboo_colors[2])
# plt.ylim([-6, 1])
# plt.xlim([-6, 1])
plt.savefig('%s/%s/normalization/%s_n-pos_vs_n-neg_ratio.pdf' % (output_dir, folder, folder))
plt.show()






