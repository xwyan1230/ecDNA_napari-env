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

folder = '5uM_48hr'
batches = [1, 2, 3]
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
fov_area = 2.336  # mm2
well_area = 11
test1 = ['XY0%s' % (x+1) for x in range(9)] + ['XY%s' % (x+10) for x in range(50)]
test2 = ['XY0%s' % (x+1) for x in range(9)] + ['XY%s' % (x+10) for x in range(60)]
ylim_val = [-500, 17500]

data = pd.DataFrame()
for batch in batches:
    df = pd.read_csv('%s/%s/%s_%s.txt' % (data_dir, folder, folder, batch), na_values=['.'], sep='\t')
    df['log10_fov_hoechst'] = np.log10(df['fov_hoechst'])
    data = pd.concat([data, df], axis=0)

data['log10_fov_hoechst_per_well'] = data['log10_fov_hoechst']*well_area/fov_area
data['n_filtered_per_well'] = data['n_filtered']*well_area/fov_area

data1 = data[(data['rep'] != 0) & (data['plate'] == batches[0]) & (data['sample'].isin(test1))].copy().reset_index(drop=True)
data2 = data[(data['rep'] != 0) & (data['plate'] == batches[0]) & (~data['sample'].isin(test1))].copy().reset_index(drop=True)
data3 = data[(data['rep'] != 0) & (data['plate'] == batches[1]) & (data['sample'].isin(test2))].copy().reset_index(drop=True)
data4 = data[(data['rep'] != 0) & (data['plate'] == batches[1]) & (~data['sample'].isin(test2))].copy().reset_index(drop=True)
data5 = data[(data['rep'] != 0) & (data['plate'] == batches[2])].copy().reset_index(drop=True)
data1_ctrl = data1[data1['treatment']=='DMSO'].copy().reset_index(drop=True)
data2_ctrl = data2[data2['treatment']=='DMSO'].copy().reset_index(drop=True)
data3_ctrl = data3[data3['treatment']=='DMSO'].copy().reset_index(drop=True)
data4_ctrl = data4[data4['treatment']=='DMSO'].copy().reset_index(drop=True)
data5_ctrl = data5[data5['treatment']=='DMSO'].copy().reset_index(drop=True)

plt.subplots(figsize=(9, 7))
feature = 'log10_fov_hoechst_per_well'
feature1 = 'n_filtered_per_well'
sns.scatterplot(data=data1, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[0])
sns.scatterplot(data=data2, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[2])
sns.scatterplot(data=data3, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[4])
sns.scatterplot(data=data4, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[6])
sns.scatterplot(data=data5, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[8])
plt.ylim(ylim_val)
plt.xlim([40, 50])
if not os.path.exists('%s/%s/normalization/' % (output_dir, folder)):
    os.makedirs('%s/%s/normalization/' % (output_dir, folder))
plt.savefig('%s/%s/normalization/%s_before.pdf' % (output_dir, folder, folder))
plt.show()

plt.subplots(figsize=(9, 7))
feature = 'log10_fov_hoechst_per_well'
feature1 = 'n_filtered_per_well'
sns.scatterplot(data=data1_ctrl, x=feature, y=feature1, alpha=1, s=30, color= rainboo_colors[0])
sns.scatterplot(data=data2_ctrl, x=feature, y=feature1, alpha=1, s=30, color= rainboo_colors[2])
sns.scatterplot(data=data3_ctrl, x=feature, y=feature1, alpha=1, s=30, color= rainboo_colors[4])
sns.scatterplot(data=data4_ctrl, x=feature, y=feature1, alpha=1, s=30, color= rainboo_colors[6])
sns.scatterplot(data=data5_ctrl, x=feature, y=feature1, alpha=1, s=30, color= rainboo_colors[8])
plt.ylim(ylim_val)
plt.xlim([40, 50])
plt.show()


def get_normalization(data, ctrl, ref):
    mean_hoechst = np.mean(ref[ref['n_filtered']>1000]['log10_fov_hoechst'])
    mean_n = np.mean(ref[ref['n_filtered']>1000]['n_filtered'])
    hoechst_ratio = np.mean(ctrl[ctrl['n_filtered']>1000]['log10_fov_hoechst'])/mean_hoechst
    n_ratio = np.mean(ctrl[ctrl['n_filtered']>1000]['n_filtered'])/mean_n
    data['log10_fov_hoechst_normalized'] = data['log10_fov_hoechst']/hoechst_ratio
    data['n_filtered_normalized'] = data['n_filtered']/n_ratio
    data['log10_fov_hoechst_normalized_per_well'] = data['log10_fov_hoechst_normalized'] * well_area / fov_area
    data['n_filtered_normalized_per_well'] = data['n_filtered_normalized'] * well_area / fov_area
    return data


data1 = get_normalization(data1, data1_ctrl, data1_ctrl)
data2 = get_normalization(data2, data2_ctrl, data1_ctrl)
data3 = get_normalization(data3, data3_ctrl, data1_ctrl)
data4 = get_normalization(data4, data4_ctrl, data1_ctrl)
data5 = get_normalization(data5, data5_ctrl, data1_ctrl)
data1_ctrl = data1[data1['treatment']=='DMSO'].copy().reset_index(drop=True)
data2_ctrl = data2[data2['treatment']=='DMSO'].copy().reset_index(drop=True)
data3_ctrl = data3[data3['treatment']=='DMSO'].copy().reset_index(drop=True)
data4_ctrl = data4[data4['treatment']=='DMSO'].copy().reset_index(drop=True)
data5_ctrl = data5[data5['treatment']=='DMSO'].copy().reset_index(drop=True)

data_all = pd.concat([data1, data2, data3, data4, data5], axis=0)
data_all.to_csv('%s/%s/%s_normalized.txt' % (data_dir, folder, folder), index=False, sep='\t')

plt.subplots(figsize=(9, 7))
feature = 'log10_fov_hoechst_normalized_per_well'
feature1 = 'n_filtered_normalized_per_well'
sns.scatterplot(data=data1, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[0])
sns.scatterplot(data=data2, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[2])
sns.scatterplot(data=data3, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[4])
sns.scatterplot(data=data4, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[6])
sns.scatterplot(data=data5, x=feature, y=feature1, alpha=1, s=20, color= rainboo_colors[8])
plt.ylim(ylim_val)
plt.xlim([40, 50])
plt.savefig('%s/%s/normalization/%s_after.pdf' % (output_dir, folder, folder))
plt.show()

plt.subplots(figsize=(9, 7))
feature = 'log10_fov_hoechst_normalized_per_well'
feature1 = 'n_filtered_normalized_per_well'
sns.scatterplot(data=data_all, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
sns.scatterplot(data=data_all[data_all['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[0])
plt.ylim([-500, 14000])
plt.xlim([40, 50])
plt.savefig('%s/%s/normalization/%s_after_DMSO.pdf' % (output_dir, folder, folder))
plt.show()

plt.subplots(figsize=(9, 7))
feature = 'log10_fov_hoechst_normalized_per_well'
feature1 = 'n_filtered_normalized_per_well'
sns.scatterplot(data=data1_ctrl, x=feature, y=feature1, alpha=1, s=30, color= rainboo_colors[0])
sns.scatterplot(data=data2_ctrl, x=feature, y=feature1, alpha=1, s=30, color= rainboo_colors[2])
sns.scatterplot(data=data3_ctrl, x=feature, y=feature1, alpha=1, s=30, color= rainboo_colors[4])
sns.scatterplot(data=data4_ctrl, x=feature, y=feature1, alpha=1, s=30, color= rainboo_colors[6])
sns.scatterplot(data=data5_ctrl, x=feature, y=feature1, alpha=1, s=30, color= rainboo_colors[8])
plt.ylim(ylim_val)
plt.xlim([40, 50])
plt.show()

print(data_all[(data_all[feature]>47)&(data_all[feature1]<2000)]['treatment'])

df = pd.DataFrame()
data_rep1 = data_all[data_all['rep']==1].copy().reset_index(drop=True).sort_values(by='group').reset_index(drop=True)
data_rep2 = data_all[data_all['rep']==2].copy().reset_index(drop=True).sort_values(by='group').reset_index(drop=True)
if data_rep1['group'].tolist() == data_rep2['group'].tolist():
    df['treatment'] = data_rep1['treatment']
    df['hoechst_rep1'] = data_rep1['log10_fov_hoechst_normalized_per_well']
    df['hoechst_rep2'] = data_rep2['log10_fov_hoechst_normalized_per_well']
    df['n_rep1'] = data_rep1['n_filtered_normalized_per_well']
    df['n_rep2'] = data_rep2['n_filtered_normalized_per_well']
    df['mean_hoechst'] = (df['hoechst_rep1']+df['hoechst_rep2'])/2
    df['mean_n'] = (df['n_rep1']+df['n_rep2'])/2

feature = 'log10_fov_hoechst_normalized_per_well'
rep1_std = np.std(data_all[(data_all['rep'] == 1) & (data_all['treatment'] == 'DMSO')][feature])
rep2_std = np.std(data_all[(data_all['rep'] == 2) & (data_all['treatment'] == 'DMSO')][feature])
rep1_mean = np.mean(data_all[(data_all['rep'] == 1) & (data_all['treatment'] == 'DMSO')][feature])
rep2_mean = np.mean(data_all[(data_all['rep'] == 2) & (data_all['treatment'] == 'DMSO')][feature])
print(rep1_std)
print(rep1_mean)
print(rep2_std)
print(rep2_mean)
df['log2_hoechst_rep1'] = np.log2(df['hoechst_rep1']/rep1_mean)
df['log2_hoechst_rep2'] = np.log2(df['hoechst_rep2']/rep2_mean)
df['log2_mean_hoechst'] = (df['log2_hoechst_rep1']+df['log2_hoechst_rep2'])/2

plt.subplots(figsize=(9, 7))
feature = 'hoechst_rep1'
feature1 = 'hoechst_rep2'
data_hoechst_hit = df[(df[feature] < (rep1_mean-3*rep1_std)) & (df[feature1] < (rep2_mean-3*rep2_std))].copy().reset_index(drop=True)
data_hoechst_hit.to_csv('%s/%s/%s_hoechst_hit.txt' % (data_dir, folder, folder), index=False, sep='\t')
sns.scatterplot(data=df, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
sns.scatterplot(data=df[df['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
sns.scatterplot(data=data_hoechst_hit, x=feature, y=feature1, alpha=1, s=20, color='#93b777')
plt.axvline(x=rep1_mean-3*rep1_std, linestyle='--', color=rainboo_colors[0])
plt.axhline(y=rep2_mean-3*rep2_std, linestyle='--', color=rainboo_colors[0])
plt.ylim([40, 50])
plt.xlim([40, 50])
plt.savefig('%s/%s/normalization/%s_hoechst_rep.pdf' % (output_dir, folder, folder))
plt.show()

plt.subplots(figsize=(9, 7))
feature = 'log2_hoechst_rep1'
feature1 = 'log2_hoechst_rep2'
# data_hoechst_hit = df[(df[feature] < (rep1_mean-3*rep1_std)) & (df[feature1] < (rep2_mean-3*rep2_std))].copy().reset_index(drop=True)
# data_hoechst_hit.to_csv('%s/%s/%s_hoechst_hit.txt' % (data_dir, folder, folder), index=False, sep='\t')
sns.scatterplot(data=df, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
sns.scatterplot(data=df[df['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
# sns.scatterplot(data=data_hoechst_hit, x=feature, y=feature1, alpha=1, s=20, color='#93b777')
# plt.axvline(x=rep1_mean-3*rep1_std, linestyle='--', color=rainboo_colors[0])
# plt.axhline(y=rep2_mean-3*rep2_std, linestyle='--', color=rainboo_colors[0])
# plt.ylim([40, 50])
# plt.xlim([40, 50])
plt.savefig('%s/%s/normalization/%s_hoechst_rep_log2.pdf' % (output_dir, folder, folder))
plt.show()

feature = 'n_filtered_normalized_per_well'
rep1_std = np.std(data_all[(data_all['rep'] == 1) & (data_all['treatment'] == 'DMSO')][feature])
rep2_std = np.std(data_all[(data_all['rep'] == 2) & (data_all['treatment'] == 'DMSO')][feature])
rep1_mean = np.mean(data_all[(data_all['rep'] == 1) & (data_all['treatment'] == 'DMSO')][feature])
rep2_mean = np.mean(data_all[(data_all['rep'] == 2) & (data_all['treatment'] == 'DMSO')][feature])
print(rep1_std)
print(rep1_mean)
print(rep2_std)
print(rep2_mean)
df['log2_n_rep1'] = np.log2(df['n_rep1']/rep1_mean)
df['log2_n_rep2'] = np.log2(df['n_rep2']/rep2_mean)
df['log2_mean_n'] = (df['log2_n_rep1']+df['log2_n_rep2'])/2
df.to_csv('%s/%s/%s_average.txt' % (data_dir, folder, folder), index=False, sep='\t')

plt.subplots(figsize=(9, 7))
feature = 'n_rep1'
feature1 = 'n_rep2'
data_n_hit = df[(df[feature] < (rep1_mean-3*rep1_std)) & (df[feature1] < (rep2_mean-3*rep2_std))].copy().reset_index(drop=True)
data_n_hit.to_csv('%s/%s/%s_n_hit.txt' % (data_dir, folder, folder), index=False, sep='\t')
sns.scatterplot(data=df, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
sns.scatterplot(data=df[df['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
sns.scatterplot(data=data_n_hit, x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[2])
plt.axvline(x=rep1_mean-3*rep1_std, linestyle='--', color=rainboo_colors[0])
plt.axhline(y=rep2_mean-3*rep2_std, linestyle='--', color=rainboo_colors[0])
plt.ylim([-500, 16000])
plt.xlim([-500, 16000])
plt.savefig('%s/%s/normalization/%s_n_rep.pdf' % (output_dir, folder, folder))
plt.show()

plt.subplots(figsize=(9, 7))
feature = 'log2_n_rep1'
feature1 = 'log2_n_rep2'
# data_hoechst_hit = df[(df[feature] < (rep1_mean-3*rep1_std)) & (df[feature1] < (rep2_mean-3*rep2_std))].copy().reset_index(drop=True)
# data_hoechst_hit.to_csv('%s/%s/%s_hoechst_hit.txt' % (data_dir, folder, folder), index=False, sep='\t')
sns.scatterplot(data=df, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
sns.scatterplot(data=df[df['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
# sns.scatterplot(data=data_hoechst_hit, x=feature, y=feature1, alpha=1, s=20, color='#93b777')
# plt.axvline(x=rep1_mean-3*rep1_std, linestyle='--', color=rainboo_colors[0])
# plt.axhline(y=rep2_mean-3*rep2_std, linestyle='--', color=rainboo_colors[0])
# plt.ylim([40, 50])
# plt.xlim([40, 50])
plt.savefig('%s/%s/normalization/%s_n_rep_log2.pdf' % (output_dir, folder, folder))
plt.show()

setA = data_hoechst_hit['treatment']
setB = data_n_hit['treatment']

print(set(setA) - set(setB))
print(set(setB) - set(setA))

plt.subplots(figsize=(9, 9))
venn2([set(setA), set(setB)], set_labels=['hoechst', 'n'], set_colors=(rainboo_colors[1], rainboo_colors[6]))
plt.savefig('%s/%s/normalization/%s_venn2_hoechst_vs_n.pdf' % (output_dir, folder, folder))
plt.show()

treatments_lst = list(set(setA) - set(setB))
treatments_lst1 = list(set(setB) - set(setA))
treatments_lst2 = list(set(setA) & set(setB))

plt.subplots(figsize=(9, 7))
feature = 'log10_fov_hoechst_normalized_per_well'
feature1 = 'n_filtered_normalized_per_well'
sns.scatterplot(data=data_all, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
sns.scatterplot(data=data_all[data_all['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
sns.scatterplot(data=data_all[data_all['treatment'].isin(treatments_lst1)], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[2])
sns.scatterplot(data=data_all[data_all['treatment'].isin(treatments_lst)], x=feature, y=feature1, alpha=1, s=20, color='#93b777')
sns.scatterplot(data=data_all[data_all['treatment'].isin(treatments_lst2)], x=feature, y=feature1, alpha=1, s=20, color='#edb08e')  # '#d84480' '#8f44d8'
plt.ylim([-500, 16000])
plt.xlim([40, 50])
plt.savefig('%s/%s/normalization/%s_additional_hits.pdf' % (output_dir, folder, folder))
plt.show()

plt.subplots(figsize=(9, 7))
feature = 'mean_hoechst'
feature1 = 'mean_n'
sns.scatterplot(data=df, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
sns.scatterplot(data=df[df['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
sns.scatterplot(data=df[df['treatment'].isin(treatments_lst1)], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[2])
sns.scatterplot(data=df[df['treatment'].isin(treatments_lst)], x=feature, y=feature1, alpha=1, s=20, color='#93b777')
sns.scatterplot(data=df[df['treatment'].isin(treatments_lst2)], x=feature, y=feature1, alpha=1, s=20, color='#edb08e')  # '#d84480' '#8f44d8'
plt.ylim([-500, 16000])
plt.xlim([40, 50])
plt.savefig('%s/%s/normalization/%s_additional_hits_average.pdf' % (output_dir, folder, folder))
plt.show()

plt.subplots(figsize=(9, 7))
feature = 'mean_hoechst'
feature1 = 'mean_n'
sns.scatterplot(data=df, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
sns.scatterplot(data=df[df['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
data_flt = df[df['treatment'].isin(treatments_lst1)].copy().reset_index(drop=True)
sns.scatterplot(data=df[df['treatment'].isin(treatments_lst1)], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[2])
sns.scatterplot(data=df[df['treatment'].isin(treatments_lst)], x=feature, y=feature1, alpha=1, s=20, color='#93b777')
sns.scatterplot(data=df[df['treatment'].isin(treatments_lst2)], x=feature, y=feature1, alpha=1, s=20, color='#edb08e')  # '#d84480' '#8f44d8'
for i in range(len(data_flt)):
    plt.text(x=data_flt[feature][i] + 0.1, y=data_flt[feature1][i], s=data_flt['treatment'][i],
             size=6, color=rainboo_colors[2])
plt.ylim([-500, 16000])
plt.xlim([40, 50])
plt.savefig('%s/%s/normalization/%s_additional_hits_average_txt.pdf' % (output_dir, folder, folder))
plt.show()



