import numpy as np
import pandas as pd
import utilities as uti
import os

# INPUT PARAMETERS
# file info
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
batch = '5uM_48hr'
total_group = 290

# PARAMETERS
data = pd.read_csv('%s/01_summary/%s_grouping.txt' % (data_dir, batch), na_values=['.'], sep='\t')
print(len(data))

feature_lst = ['batch', 'group', 'treatment', 'cc_filter', 'qc_filter', 'cytokinesis_hit']
feature_lst1 = ['log2fc_hoechst_mean']
feature_lst2 = ['log2fc_n_mean', 'log2fc_n_neg_mean', 'log2fc_n_pos_mean', 'log2fc_ratio_mean']
feature_lst3 = ['per_G1_mean', 'per_G1S_mean', 'per_S_mean', 'per_G2M_mean', 'per_G2MG1_mean',
                           'per_neg_G1_mean', 'per_neg_G1S_mean', 'per_neg_S_mean', 'per_neg_G2M_mean',
                           'per_neg_G2MG1_mean',
                           'per_pos_G1_mean', 'per_pos_G1S_mean', 'per_pos_S_mean', 'per_pos_G2M_mean',
                           'per_pos_G2MG1_mean',
                           'hoechst_G1_mean', 'hoechst_G1S_mean', 'hoechst_S_mean', 'hoechst_G2M_mean', 'hoechst_G2MG1_mean',
                           'hoechst_neg_G1_mean', 'hoechst_neg_G1S_mean', 'hoechst_neg_S_mean', 'hoechst_neg_G2M_mean',
                           'hoechst_neg_G2MG1_mean',
                           'hoechst_pos_G1_mean', 'hoechst_pos_G1S_mean', 'hoechst_pos_S_mean', 'hoechst_pos_G2M_mean',
                           'hoechst_pos_G2MG1_mean']
df = pd.DataFrame()
for feature in feature_lst:
    print(feature)
    df[feature] = data[feature]
for feature in feature_lst1:
    print(feature)
    df[feature] = uti.call_hit(data, feature, ['qc_filter'], [])
for feature in feature_lst2:
    print(feature)
    df[feature] = uti.call_hit(data, feature, ['qc_filter'], ['cytokinesis_hit'])
for feature in feature_lst3:
    print(feature)
    df[feature] = uti.call_hit(data, feature, ['qc_filter', 'cc_filter'], [])

df.to_csv('%s/01_summary/%s_hitcalling.txt' % (output_dir, batch), index=False, sep='\t')
print("DONE!")