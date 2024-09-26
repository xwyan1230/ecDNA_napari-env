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
colors = ['#bc4d4a', '#669daa', '#dd933c'] * 3
fov_area = 2.336  # mm2
well_area = 11
ylim_val = [-500, 17500]

data = pd.DataFrame()
for batch in batches:
    df = pd.read_csv('%s/%s/%s_%s_update.txt' % (data_dir, folder, folder, batch), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0)

features = ['n_G1', 'n_S', 'n_G2M']
feature_sum = ['n_filtered'] * 3
for i in range(len(features)):
    data['per_%s' % features[i]] = data[features[i]]/data[feature_sum[i]]


def model_func(x, k, b):
    return np.exp(k*np.log(x)+b)


def cc_normalize(data, phase):
    if phase == 'G1':
        mean_val = 0.34750426368895027
        k = -0.5101334598046143
        b = -0.7021534884918927
    elif phase == 'S':
        mean_val = 0.30149297837435207
        k = -0.5112208788698921
        b = -0.7330796350779476
    else:
        mean_val = 0.1253451744762164
        k = -0.48540086608779354
        b = -1.2088281578713023
    std_val = model_func(data['n_filtered'], k, b)
    data['per_%s_nor' % phase] = (mean_val-data['per_n_%s' % phase])/-std_val
    return data


data = cc_normalize(data, 'G1')
data = cc_normalize(data, 'S')
data = cc_normalize(data, 'G2M')


def get_average(data, data_rep1, data_rep2, feature):
    data['%s_rep1' % feature] = data_rep1[feature]
    data['%s_rep2' % feature] = data_rep2[feature]
    data['mean_%s' % feature] = (data['%s_rep1' % feature] + data['%s_rep2' % feature])/2
    return data


df = pd.DataFrame()
data_rep1 = data[data['rep']==1].copy().reset_index(drop=True).sort_values(by='group').reset_index(drop=True)
data_rep2 = data[data['rep']==2].copy().reset_index(drop=True).sort_values(by='group').reset_index(drop=True)
features = ['per_G1_nor', 'per_S_nor', 'per_G2M_nor']
if data_rep1['group'].tolist() == data_rep2['group'].tolist():
    df['treatment'] = data_rep1['treatment']
    for feature in features:
        df = get_average(df, data_rep1, data_rep2, feature)
df.to_csv('%s/%s/%s_cellcycle_average_normalize.txt' % (data_dir, folder, folder), index=False, sep='\t')

for k in range(len(features)):
    f = features[k]
    rep1_std = np.std(data[(data['rep'] == 1) & (data['treatment'] == 'DMSO')][f])
    rep2_std = np.std(data[(data['rep'] == 2) & (data['treatment'] == 'DMSO')][f])
    rep1_mean = np.mean(data[(data['rep'] == 1) & (data['treatment'] == 'DMSO')][f])
    rep2_mean = np.mean(data[(data['rep'] == 2) & (data['treatment'] == 'DMSO')][f])
    print(rep1_std)
    print(rep1_mean)
    print(rep2_std)
    print(rep2_mean)

    plt.subplots(figsize=(9, 7))
    feature = '%s_rep1' % f
    feature1 = '%s_rep2' % f
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
    plt.savefig('%s/%s/cellcycle/%s_%s_rep_normalize.pdf' % (output_dir, folder, folder, feature))
    plt.show()








