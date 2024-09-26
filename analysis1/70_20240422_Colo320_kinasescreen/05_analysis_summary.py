import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import seaborn as sns
from scipy.stats import ks_2samp, ttest_ind
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sprocessed/" % master_folder

batch = '5uM_24hr'
plates = [1, 2, 3]
rep = [1,2]
cell = 'DF+HFH'

target_df = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % master_folder, na_values=['.'])
print(target_df.head())

data = pd.DataFrame()
for plate in plates:
    df = pd.read_csv('%s/%s/%s_%s.txt' % (data_dir, batch, batch, plate), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0)
data_screen = data[(data['cell'] == cell) & (data['rep'].isin(rep))].copy().reset_index(drop=True)
data_screen['pos_vs_neg'] = (data_screen['n_pos']+0.01)/(data_screen['n_neg']+0.01)
data_screen['cytokinesis'] = data_screen['fov_hoechst']/data_screen['n_filtered']

print(max(data_screen['group']))

data_summary = pd.DataFrame(columns=['screen', 'group', 'cell', 'treatment', 'target',
                                     'mean_n_filtered', 'log2fc_n_filtered', 'p_n_filtered', 'mlog10p_n_filtered',
                                     'mean_n_neg', 'log2fc_n_neg', 'p_n_neg', 'mlog10p_n_neg',
                                     'mean_n_pos', 'log2fc_n_pos', 'p_n_pos', 'mlog10p_n_pos',
                                     'mean_fov_hoechst', 'log2fc_fov_hoechst', 'p_fov_hoechst', 'mlog10p_fov_hoechst',
                                     'mean_pos_vs_neg', 'log2fc_pos_vs_neg', 'p_pos_vs_neg', 'mlog10p_pos_vs_neg',
                                     'mean_cytokinesis', 'log2fc_cytokinesis', 'p_cytokinesis', 'mlog10p_cytokinesis'])

data_ctrl = data_screen[data_screen['treatment'] == 'DMSO'].copy().reset_index(drop=True)

mean_feature_ctrls = []
features = ['n_filtered', 'n_neg', 'n_pos', 'fov_hoechst', 'pos_vs_neg', 'cytokinesis']
for feature in features:
    mean_feature_ctrls = mean_feature_ctrls + [np.mean(data_ctrl[feature])]


def get_score(data_temp, feature, mean_feature_ctrl, data_ctrl):
    mean_feature = np.mean(data_temp[feature])
    log2fc_feature = np.log2(mean_feature/mean_feature_ctrl)
    p_feature = ttest_ind(data_temp[feature].tolist(), data_ctrl[feature].tolist())[1]
    mlog10p_feature = -np.log10(p_feature)
    return mean_feature, log2fc_feature, p_feature, mlog10p_feature


for i in range(max(data_screen['group'])):
    data_temp = data_screen[data_screen['group'] == i+1].copy().reset_index(drop=True)
    treatment = data_temp['treatment'][0]
    if treatment == 'DMSO':
        target = 'DMSO'
    else:
        target = target_df[target_df['compound name'] == treatment]['TargetGene'].tolist()[0]
    feature_lst = []
    for j in range(len(features)):
        feature = features[j]
        mean_feature_ctrl = mean_feature_ctrls[j]
        mean_feature, log2fc_feature, p_feature, mlog10p_feature = get_score(data_temp, feature, mean_feature_ctrl, data_ctrl)
        feature_lst = feature_lst + [mean_feature, log2fc_feature, p_feature, mlog10p_feature]
    data_summary.loc[len(data_summary.index)] = [data_temp['screen'][0], data_temp['group'][0], data_temp['cell'][0], treatment, target] + feature_lst

data_summary.to_csv('%s/%s/%s_summary.txt' % (output_dir, batch, batch), index=False, sep='\t')

