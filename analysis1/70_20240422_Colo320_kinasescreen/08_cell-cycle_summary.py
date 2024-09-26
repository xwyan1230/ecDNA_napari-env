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

batch = '5uM_48hr'
plates = [1, 2, 3]
rep = [1,2]
cell = 'DF+HFH'

target_df = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % master_folder, na_values=['.'])
print(target_df.head())

data = pd.DataFrame()
for plate in plates:
    df = pd.read_csv('%s/%s/%s_%s_cc.txt' % (data_dir, batch, batch, plate), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0)
data_screen = data[(data['cell'] == cell) & (data['rep'].isin(rep))].copy().reset_index(drop=True)
data_screen['pos_vs_neg'] = (data_screen['n_pos']+0.01)/(data_screen['n_neg']+0.01)

print(max(data_screen['group']))

data_summary = pd.DataFrame(columns=['screen', 'group', 'cell', 'treatment', 'target',
                                     'rep', 'cc_score',
                                     'rep_neg', 'rep_pos',
                                     'cc_score_neg', 'cc_score_pos',
                                     'mean_per_G1', 'log2fc_per_G1', 'p_per_G1', 'mlog10p_per_G1',
                                     'mean_per_G1S', 'log2fc_per_G1S', 'p_per_G1S', 'mlog10p_per_G1S',
                                     'mean_per_S', 'log2fc_per_S', 'p_per_S', 'mlog10p_per_S',
                                     'mean_per_G2M', 'log2fc_per_G2M', 'p_per_G2M', 'mlog10p_per_G2M',
                                     'mean_per_neg_G1', 'log2fc_per_neg_G1', 'p_per_neg_G1', 'mlog10p_per_neg_G1',
                                     'mean_per_neg_G1S', 'log2fc_per_neg_G1S', 'p_per_neg_G1S', 'mlog10p_per_neg_G1S',
                                     'mean_per_neg_S', 'log2fc_per_neg_S', 'p_per_neg_S', 'mlog10p_per_neg_S',
                                     'mean_per_neg_G2M', 'log2fc_per_neg_G2M', 'p_per_neg_G2M', 'mlog10p_per_neg_G2M',
                                     'mean_per_pos_G1', 'log2fc_per_pos_G1', 'p_per_pos_G1', 'mlog10p_per_pos_G1',
                                     'mean_per_pos_G1S', 'log2fc_per_pos_G1S', 'p_per_pos_G1S', 'mlog10p_per_pos_G1S',
                                     'mean_per_pos_S', 'log2fc_per_pos_S', 'p_per_pos_S', 'mlog10p_per_pos_S',
                                     'mean_per_pos_G2M', 'log2fc_per_pos_G2M', 'p_per_pos_G2M', 'mlog10p_per_pos_G2M'
                                     ])

data_ctrl = data_screen[data_screen['treatment'] == 'DMSO'].copy().reset_index(drop=True)
mean_feature_ctrls = []
features = ['per_G1', 'per_G1S', 'per_S', 'per_G2M', 'per_neg_G1', 'per_neg_G1S', 'per_neg_S', 'per_neg_G2M', 'per_pos_G1', 'per_pos_G1S', 'per_pos_S', 'per_pos_G2M']
for feature in features:
    mean_feature_ctrls = mean_feature_ctrls + [np.mean(data_ctrl[feature])]


def get_score(data_temp, feature, mean_feature_ctrl, data_ctrl):
    mean_feature = np.mean(data_temp[feature])
    if mean_feature != 0:
        log2fc_feature = np.log2(mean_feature/mean_feature_ctrl)
    else:
        log2fc_feature = np.log2(0.00001 / mean_feature_ctrl)
    p_feature = ttest_ind(data_temp[feature].tolist(), data_ctrl[feature].tolist())[1]
    mlog10p_feature = -np.log10(p_feature)
    return mean_feature, log2fc_feature, p_feature, mlog10p_feature


for i in range(max(data_screen['group'])):
    data_temp = data_screen[data_screen['group'] == i+1].copy().reset_index(drop=True)
    data_temp_neg = data_temp[data_temp['n_neg'] != 0].copy().reset_index(drop=True)
    data_temp_pos = data_temp[data_temp['n_pos'] != 0].copy().reset_index(drop=True)
    treatment = data_temp['treatment'][0]
    if treatment == 'DMSO':
        target = 'DMSO'
    else:
        target = target_df[target_df['compound name'] == treatment]['TargetGene'].tolist()[0]
    feature_lst = []
    rep = len(data_temp)
    rep_neg = len(data_temp_neg)
    rep_pos = len(data_temp_pos)
    for j in range(len(features)):
        feature = features[j]
        mean_feature_ctrl = mean_feature_ctrls[j]
        if feature[4] == 'n':
            mean_feature, log2fc_feature, p_feature, mlog10p_feature = get_score(data_temp_neg, feature, mean_feature_ctrl, data_ctrl)
        elif feature[4] == 'p':
            mean_feature, log2fc_feature, p_feature, mlog10p_feature = get_score(data_temp_pos, feature,
                                                                                 mean_feature_ctrl, data_ctrl)
        else:
            mean_feature, log2fc_feature, p_feature, mlog10p_feature = get_score(data_temp, feature,
                                                                                 mean_feature_ctrl, data_ctrl)
        feature_lst = feature_lst + [mean_feature, log2fc_feature, p_feature, mlog10p_feature]
    cc_score = np.abs(feature_lst[1]) + np.abs(feature_lst[5]) + np.abs(feature_lst[9]) + np.abs(feature_lst[13])
    cc_score_neg = np.abs(feature_lst[17]) + np.abs(feature_lst[21]) + np.abs(feature_lst[25]) + np.abs(feature_lst[29])
    cc_score_pos = np.abs(feature_lst[33]) + np.abs(feature_lst[37]) + np.abs(feature_lst[41]) + np.abs(feature_lst[45])

    data_summary.loc[len(data_summary.index)] = [data_temp['screen'][0], data_temp['group'][0], data_temp['cell'][0],
                                                 treatment, target, rep, cc_score, rep_neg, rep_pos, cc_score_neg, cc_score_pos] + feature_lst

mean_feature_ctrls = []
data_summary_ctrl = data_summary[data_summary['treatment'] == 'DMSO'].copy().reset_index(drop=True)
features = ['cc_score', 'cc_score_neg', 'cc_score_pos']
for feature in features:
    mean_feature_ctrls = mean_feature_ctrls + [np.mean(data_summary_ctrl[feature])]
print(mean_feature_ctrls)
data_summary['cc_score'] = data_summary['cc_score'] - mean_feature_ctrls[0]
data_summary['cc_score_neg'] = data_summary['cc_score_neg'] - mean_feature_ctrls[1]
data_summary['cc_score_pos'] = data_summary['cc_score_pos'] - mean_feature_ctrls[2]

data_summary.to_csv('%s/%s/%s_summary_cc.txt' % (output_dir, batch, batch), index=False, sep='\t')

