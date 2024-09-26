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
data = pd.read_csv('%s/01_summary/%s_normalize.txt' % (data_dir, batch), na_values=['.'], sep='\t')
print(len(data))
data['ratio_HSR_DM'] = data['n_pos']/(data['n_neg']+0.1)
data = uti.get_cc_percentage(data)
data_cc = data[data['cc_filter'] == 1].copy().reset_index(drop=True)
print(len(data_cc))

# hoechst: log10_fov_hoechst_normalized
# n: n_filtered_normalized
# n_neg: n_neg_normalized
# n_pos: n_pos_normalized
features = ['log10_fov_hoechst_normalized', 'n_filtered_normalized', 'n_neg_normalized', 'n_pos_normalized']

df = pd.DataFrame(columns=['batch', 'group', 'treatment',
                           'hoechst_rep1', 'hoechst_rep2', 'hoechst_mean',
                           'log2fc_hoechst_rep1', 'log2fc_hoechst_rep2', 'log2fc_hoechst_mean',
                           'n_rep1', 'n_rep2', 'n_mean',
                           'log2fc_n_rep1', 'log2fc_n_rep2', 'log2fc_n_mean',
                           'n_neg_rep1', 'n_neg_rep2', 'n_neg_mean',
                           'log2fc_n_neg_rep1', 'log2fc_n_neg_rep2', 'log2fc_n_neg_mean',
                           'n_pos_rep1', 'n_pos_rep2', 'n_pos_mean',
                           'log2fc_n_pos_rep1', 'log2fc_n_pos_rep2', 'log2fc_n_pos_mean',
                           'ratio_rep1', 'ratio_rep2', 'ratio_mean',
                           'log2fc_ratio_rep1', 'log2fc_ratio_rep2', 'log2fc_ratio_mean',
                           'per_G1_rep1', 'per_G1_rep2', 'per_G1_mean',
                           'per_G1S_rep1', 'per_G1S_rep2', 'per_G1S_mean',
                           'per_S_rep1', 'per_S_rep2', 'per_S_mean',
                           'per_G2M_rep1', 'per_G2M_rep2', 'per_G2M_mean',
                           'per_G2MG1_rep1', 'per_G2MG1_rep2', 'per_G2MG1_mean',
                           'delta_G1_rep1', 'delta_G1_rep2', 'delta_G1_mean',
                           'delta_G1S_rep1', 'delta_G1S_rep2', 'delta_G1S_mean',
                           'delta_S_rep1', 'delta_S_rep2', 'delta_S_mean',
                           'delta_G2M_rep1', 'delta_G2M_rep2', 'delta_G2M_mean',
                           'delta_G2MG1_rep1', 'delta_G2MG1_rep2', 'delta_G2MG1_mean',
                           'per_neg_G1_rep1', 'per_neg_G1_rep2', 'per_neg_G1_mean',
                           'per_neg_G1S_rep1', 'per_neg_G1S_rep2', 'per_neg_G1S_mean',
                           'per_neg_S_rep1', 'per_neg_S_rep2', 'per_neg_S_mean',
                           'per_neg_G2M_rep1', 'per_neg_G2M_rep2', 'per_neg_G2M_mean',
                           'per_neg_G2MG1_rep1', 'per_neg_G2MG1_rep2', 'per_neg_G2MG1_mean',
                           'delta_neg_G1_rep1', 'delta_neg_G1_rep2', 'delta_neg_G1_mean',
                           'delta_neg_G1S_rep1', 'delta_neg_G1S_rep2', 'delta_neg_G1S_mean',
                           'delta_neg_S_rep1', 'delta_neg_S_rep2', 'delta_neg_S_mean',
                           'delta_neg_G2M_rep1', 'delta_neg_G2M_rep2', 'delta_neg_G2M_mean',
                           'delta_neg_G2MG1_rep1', 'delta_neg_G2MG1_rep2', 'delta_neg_G2MG1_mean',
                           'per_pos_G1_rep1', 'per_pos_G1_rep2', 'per_pos_G1_mean',
                           'per_pos_G1S_rep1', 'per_pos_G1S_rep2', 'per_pos_G1S_mean',
                           'per_pos_S_rep1', 'per_pos_S_rep2', 'per_pos_S_mean',
                           'per_pos_G2M_rep1', 'per_pos_G2M_rep2', 'per_pos_G2M_mean',
                           'per_pos_G2MG1_rep1', 'per_pos_G2MG1_rep2', 'per_pos_G2MG1_mean',
                           'delta_pos_G1_rep1', 'delta_pos_G1_rep2', 'delta_pos_G1_mean',
                           'delta_pos_G1S_rep1', 'delta_pos_G1S_rep2', 'delta_pos_G1S_mean',
                           'delta_pos_S_rep1', 'delta_pos_S_rep2', 'delta_pos_S_mean',
                           'delta_pos_G2M_rep1', 'delta_pos_G2M_rep2', 'delta_pos_G2M_mean',
                           'delta_pos_G2MG1_rep1', 'delta_pos_G2MG1_rep2', 'delta_pos_G2MG1_mean',
                           'delta_sum_rep1', 'delta_sum_rep2', 'delta_sum_mean',
                           'delta_neg_sum_rep1', 'delta_neg_sum_rep2', 'delta_neg_sum_mean',
                           'delta_pos_sum_rep1', 'delta_pos_sum_rep2', 'delta_pos_sum_mean', 'cc_filter',
                           'hoechst_G1_rep1', 'hoechst_G1_rep2', 'hoechst_G1_mean',
                           'hoechst_G1S_rep1', 'hoechst_G1S_rep2', 'hoechst_G1S_mean',
                           'hoechst_S_rep1', 'hoechst_S_rep2', 'hoechst_S_mean',
                           'hoechst_G2M_rep1', 'hoechst_G2M_rep2', 'hoechst_G2M_mean',
                           'hoechst_G2MG1_rep1', 'hoechst_G2MG1_rep2', 'hoechst_G2MG1_mean',
                           'hoechst_neg_G1_rep1', 'hoechst_neg_G1_rep2', 'hoechst_neg_G1_mean',
                           'hoechst_neg_G1S_rep1', 'hoechst_neg_G1S_rep2', 'hoechst_neg_G1S_mean',
                           'hoechst_neg_S_rep1', 'hoechst_neg_S_rep2', 'hoechst_neg_S_mean',
                           'hoechst_neg_G2M_rep1', 'hoechst_neg_G2M_rep2', 'hoechst_neg_G2M_mean',
                           'hoechst_neg_G2MG1_rep1', 'hoechst_neg_G2MG1_rep2', 'hoechst_neg_G2MG1_mean',
                           'hoechst_pos_G1_rep1', 'hoechst_pos_G1_rep2', 'hoechst_pos_G1_mean',
                           'hoechst_pos_G1S_rep1', 'hoechst_pos_G1S_rep2', 'hoechst_pos_G1S_mean',
                           'hoechst_pos_S_rep1', 'hoechst_pos_S_rep2', 'hoechst_pos_S_mean',
                           'hoechst_pos_G2M_rep1', 'hoechst_pos_G2M_rep2', 'hoechst_pos_G2M_mean',
                           'hoechst_pos_G2MG1_rep1', 'hoechst_pos_G2MG1_rep2', 'hoechst_pos_G2MG1_mean',
                           'log2fc_hoechst_G1_rep1', 'log2fc_hoechst_G1_rep2', 'log2fc_hoechst_G1_mean',
                           'log2fc_hoechst_G1S_rep1', 'log2fc_hoechst_G1S_rep2', 'log2fc_hoechst_G1S_mean',
                           'log2fc_hoechst_S_rep1', 'log2fc_hoechst_S_rep2', 'log2fc_hoechst_S_mean',
                           'log2fc_hoechst_G2M_rep1', 'log2fc_hoechst_G2M_rep2', 'log2fc_hoechst_G2M_mean',
                           'log2fc_hoechst_G2MG1_rep1', 'log2fc_hoechst_G2MG1_rep2', 'log2fc_hoechst_G2MG1_mean',
                           'log2fc_hoechst_neg_G1_rep1', 'log2fc_hoechst_neg_G1_rep2', 'log2fc_hoechst_neg_G1_mean',
                           'log2fc_hoechst_neg_G1S_rep1', 'log2fc_hoechst_neg_G1S_rep2', 'log2fc_hoechst_neg_G1S_mean',
                           'log2fc_hoechst_neg_S_rep1', 'log2fc_hoechst_neg_S_rep2', 'log2fc_hoechst_neg_S_mean',
                           'log2fc_hoechst_neg_G2M_rep1', 'log2fc_hoechst_neg_G2M_rep2', 'log2fc_hoechst_neg_G2M_mean',
                           'log2fc_hoechst_neg_G2MG1_rep1', 'log2fc_hoechst_neg_G2MG1_rep2', 'log2fc_hoechst_neg_G2MG1_mean',
                           'log2fc_hoechst_pos_G1_rep1', 'log2fc_hoechst_pos_G1_rep2', 'log2fc_hoechst_pos_G1_mean',
                           'log2fc_hoechst_pos_G1S_rep1', 'log2fc_hoechst_pos_G1S_rep2', 'log2fc_hoechst_pos_G1S_mean',
                           'log2fc_hoechst_pos_S_rep1', 'log2fc_hoechst_pos_S_rep2', 'log2fc_hoechst_pos_S_mean',
                           'log2fc_hoechst_pos_G2M_rep1', 'log2fc_hoechst_pos_G2M_rep2', 'log2fc_hoechst_pos_G2M_mean',
                           'log2fc_hoechst_pos_G2MG1_rep1', 'log2fc_hoechst_pos_G2MG1_rep2', 'log2fc_hoechst_pos_G2MG1_mean'
                           ])

for i in range(total_group):
    data_group = data[data['group'] == i+1].copy().reset_index(drop=True)
    if len(data_group) == 2:
        hoechst_rep1, hoechst_rep2, hoechst_mean, _, _, _, log2fc_hoechst_rep1, log2fc_hoechst_rep2, log2fc_hoechst_mean = uti.get_average_delta_log2fc(data_group, data, 'log10_fov_hoechst_normalized')
        n_rep1, n_rep2, n_mean, _, _, _, log2fc_n_rep1, log2fc_n_rep2, log2fc_n_mean = uti.get_average_delta_log2fc(data_group, data, 'n_filtered_normalized')
        n_neg_rep1, n_neg_rep2, n_neg_mean, _, _, _, log2fc_n_neg_rep1, log2fc_n_neg_rep2, log2fc_n_neg_mean = uti.get_average_delta_log2fc(data_group, data, 'n_neg_normalized')
        n_pos_rep1, n_pos_rep2, n_pos_mean, _, _, _, log2fc_n_pos_rep1, log2fc_n_pos_rep2, log2fc_n_pos_mean = uti.get_average_delta_log2fc(data_group, data, 'n_pos_normalized')
        ratio_rep1, ratio_rep2, ratio_mean, _, _, _, log2fc_ratio_rep1, log2fc_ratio_rep2, log2fc_ratio_mean = uti.get_average_delta_log2fc(data_group, data, 'ratio_HSR_DM')

        per_G1_rep1, per_G1_rep2, per_G1_mean, delta_G1_rep1, delta_G1_rep2, delta_G1_mean, _, _, _ = uti.get_average_delta_log2fc(data_group, data_cc, 'per_G1')
        per_G1S_rep1, per_G1S_rep2, per_G1S_mean, delta_G1S_rep1, delta_G1S_rep2, delta_G1S_mean, _, _, _ = uti.get_average_delta_log2fc(data_group, data_cc, 'per_G1S')
        per_S_rep1, per_S_rep2, per_S_mean, delta_S_rep1, delta_S_rep2, delta_S_mean, _, _, _ = uti.get_average_delta_log2fc(data_group, data_cc, 'per_S')
        per_G2M_rep1, per_G2M_rep2, per_G2M_mean, delta_G2M_rep1, delta_G2M_rep2, delta_G2M_mean, _, _, _ = uti.get_average_delta_log2fc(data_group, data_cc, 'per_G2M')
        per_G2MG1_rep1, per_G2MG1_rep2, per_G2MG1_mean, delta_G2MG1_rep1, delta_G2MG1_rep2, delta_G2MG1_mean, _, _, _ = uti.get_average_delta_log2fc(data_group, data_cc, 'per_G2MG1')
        delta_sum_rep1 = np.abs(delta_G1_rep1) + np.abs(delta_G1S_rep1) + np.abs(delta_S_rep1) + np.abs(delta_G2M_rep1) + np.abs(delta_G2MG1_rep1)
        delta_sum_rep2 = np.abs(delta_G1_rep2) + np.abs(delta_G1S_rep2) + np.abs(delta_S_rep2) + np.abs(delta_G2M_rep2) + np.abs(delta_G2MG1_rep2)
        delta_sum_mean = np.abs(delta_G1_mean) + np.abs(delta_G1S_mean) + np.abs(delta_S_mean) + np.abs(delta_G2M_mean) + np.abs(delta_G2MG1_mean)

        per_neg_G1_rep1, per_neg_G1_rep2, per_neg_G1_mean, delta_neg_G1_rep1, delta_neg_G1_rep2, delta_neg_G1_mean, _, _, _ = uti.get_average_delta_log2fc(data_group, data_cc, 'per_neg_G1')
        per_neg_G1S_rep1, per_neg_G1S_rep2, per_neg_G1S_mean, delta_neg_G1S_rep1, delta_neg_G1S_rep2, delta_neg_G1S_mean, _, _, _ = uti.get_average_delta_log2fc(data_group, data_cc, 'per_neg_G1S')
        per_neg_S_rep1, per_neg_S_rep2, per_neg_S_mean, delta_neg_S_rep1, delta_neg_S_rep2, delta_neg_S_mean, _, _, _ = uti.get_average_delta_log2fc(data_group, data_cc, 'per_neg_S')
        per_neg_G2M_rep1, per_neg_G2M_rep2, per_neg_G2M_mean, delta_neg_G2M_rep1, delta_neg_G2M_rep2, delta_neg_G2M_mean, _, _, _ = uti.get_average_delta_log2fc(data_group, data_cc, 'per_neg_G2M')
        per_neg_G2MG1_rep1, per_neg_G2MG1_rep2, per_neg_G2MG1_mean, delta_neg_G2MG1_rep1, delta_neg_G2MG1_rep2, delta_neg_G2MG1_mean, _, _, _ = uti.get_average_delta_log2fc(data_group, data_cc, 'per_neg_G2MG1')
        delta_neg_sum_rep1 = np.abs(delta_neg_G1_rep1) + np.abs(delta_neg_G1S_rep1) + np.abs(delta_neg_S_rep1) + np.abs(delta_neg_G2M_rep1) + np.abs(delta_neg_G2MG1_rep1)
        delta_neg_sum_rep2 = np.abs(delta_neg_G1_rep2) + np.abs(delta_neg_G1S_rep2) + np.abs(delta_neg_S_rep2) + np.abs(delta_neg_G2M_rep2) + np.abs(delta_neg_G2MG1_rep2)
        delta_neg_sum_mean = np.abs(delta_neg_G1_mean) + np.abs(delta_neg_G1S_mean) + np.abs(delta_neg_S_mean) + np.abs(delta_neg_G2M_mean) + np.abs(delta_neg_G2MG1_mean)

        per_pos_G1_rep1, per_pos_G1_rep2, per_pos_G1_mean, delta_pos_G1_rep1, delta_pos_G1_rep2, delta_pos_G1_mean, _, _, _ = uti.get_average_delta_log2fc(data_group, data_cc, 'per_pos_G1')
        per_pos_G1S_rep1, per_pos_G1S_rep2, per_pos_G1S_mean, delta_pos_G1S_rep1, delta_pos_G1S_rep2, delta_pos_G1S_mean, _, _, _ = uti.get_average_delta_log2fc(data_group, data_cc, 'per_pos_G1S')
        per_pos_S_rep1, per_pos_S_rep2, per_pos_S_mean, delta_pos_S_rep1, delta_pos_S_rep2, delta_pos_S_mean, _, _, _ = uti.get_average_delta_log2fc(data_group, data_cc, 'per_pos_S')
        per_pos_G2M_rep1, per_pos_G2M_rep2, per_pos_G2M_mean, delta_pos_G2M_rep1, delta_pos_G2M_rep2, delta_pos_G2M_mean, _, _, _ = uti.get_average_delta_log2fc(data_group, data_cc, 'per_pos_G2M')
        per_pos_G2MG1_rep1, per_pos_G2MG1_rep2, per_pos_G2MG1_mean, delta_pos_G2MG1_rep1, delta_pos_G2MG1_rep2, delta_pos_G2MG1_mean, _, _, _ = uti.get_average_delta_log2fc(data_group, data_cc, 'per_pos_G2MG1')
        delta_pos_sum_rep1 = np.abs(delta_pos_G1_rep1) + np.abs(delta_pos_G1S_rep1) + np.abs(delta_pos_S_rep1) + np.abs(delta_pos_G2M_rep1) + np.abs(delta_pos_G2MG1_rep1)
        delta_pos_sum_rep2 = np.abs(delta_pos_G1_rep2) + np.abs(delta_pos_G1S_rep2) + np.abs(delta_pos_S_rep2) + np.abs(delta_pos_G2M_rep2) + np.abs(delta_pos_G2MG1_rep2)
        delta_pos_sum_mean = np.abs(delta_pos_G1_mean) + np.abs(delta_pos_G1S_mean) + np.abs(delta_pos_S_mean) + np.abs(delta_pos_G2M_mean) + np.abs(delta_pos_G2MG1_mean)

        if len(data_cc[data_cc['group'] == i+1]) == 2:
            cc_filter = 1
        else:
            cc_filter = 0

        hoechst_G1_rep1, hoechst_G1_rep2, hoechst_G1_mean, _, _, _, log2fc_hoechst_G1_rep1, log2fc_hoechst_G1_rep2, log2fc_hoechst_G1_mean = uti.get_average_delta_log2fc(data_group, data_cc, 'hoechst_G1_normalized')
        hoechst_G1S_rep1, hoechst_G1S_rep2, hoechst_G1S_mean, _, _, _, log2fc_hoechst_G1S_rep1, log2fc_hoechst_G1S_rep2, log2fc_hoechst_G1S_mean = uti.get_average_delta_log2fc(data_group, data_cc, 'hoechst_G1S_normalized')
        hoechst_S_rep1, hoechst_S_rep2, hoechst_S_mean, _, _, _, log2fc_hoechst_S_rep1, log2fc_hoechst_S_rep2, log2fc_hoechst_S_mean = uti.get_average_delta_log2fc(data_group, data_cc, 'hoechst_S_normalized')
        hoechst_G2M_rep1, hoechst_G2M_rep2, hoechst_G2M_mean, _, _, _, log2fc_hoechst_G2M_rep1, log2fc_hoechst_G2M_rep2, log2fc_hoechst_G2M_mean = uti.get_average_delta_log2fc(data_group, data_cc, 'hoechst_G2M_normalized')
        hoechst_G2MG1_rep1, hoechst_G2MG1_rep2, hoechst_G2MG1_mean, _, _, _, log2fc_hoechst_G2MG1_rep1, log2fc_hoechst_G2MG1_rep2, log2fc_hoechst_G2MG1_mean = uti.get_average_delta_log2fc(data_group, data_cc, 'hoechst_G2MG1_normalized')
        hoechst_neg_G1_rep1, hoechst_neg_G1_rep2, hoechst_neg_G1_mean, _, _, _, log2fc_hoechst_neg_G1_rep1, log2fc_hoechst_neg_G1_rep2, log2fc_hoechst_neg_G1_mean = uti.get_average_delta_log2fc(data_group, data_cc, 'hoechst_neg_G1_normalized')
        hoechst_neg_G1S_rep1, hoechst_neg_G1S_rep2, hoechst_neg_G1S_mean, _, _, _, log2fc_hoechst_neg_G1S_rep1, log2fc_hoechst_neg_G1S_rep2, log2fc_hoechst_neg_G1S_mean = uti.get_average_delta_log2fc(data_group, data_cc, 'hoechst_neg_G1S_normalized')
        hoechst_neg_S_rep1, hoechst_neg_S_rep2, hoechst_neg_S_mean, _, _, _, log2fc_hoechst_neg_S_rep1, log2fc_hoechst_neg_S_rep2, log2fc_hoechst_neg_S_mean = uti.get_average_delta_log2fc(data_group, data_cc, 'hoechst_neg_S_normalized')
        hoechst_neg_G2M_rep1, hoechst_neg_G2M_rep2, hoechst_neg_G2M_mean, _, _, _, log2fc_hoechst_neg_G2M_rep1, log2fc_hoechst_neg_G2M_rep2, log2fc_hoechst_neg_G2M_mean = uti.get_average_delta_log2fc(data_group, data_cc, 'hoechst_neg_G2M_normalized')
        hoechst_neg_G2MG1_rep1, hoechst_neg_G2MG1_rep2, hoechst_neg_G2MG1_mean, _, _, _, log2fc_hoechst_neg_G2MG1_rep1, log2fc_hoechst_neg_G2MG1_rep2, log2fc_hoechst_neg_G2MG1_mean = uti.get_average_delta_log2fc(data_group, data_cc, 'hoechst_neg_G2MG1_normalized')
        hoechst_pos_G1_rep1, hoechst_pos_G1_rep2, hoechst_pos_G1_mean, _, _, _, log2fc_hoechst_pos_G1_rep1, log2fc_hoechst_pos_G1_rep2, log2fc_hoechst_pos_G1_mean = uti.get_average_delta_log2fc(data_group, data_cc, 'hoechst_pos_G1_normalized')
        hoechst_pos_G1S_rep1, hoechst_pos_G1S_rep2, hoechst_pos_G1S_mean, _, _, _, log2fc_hoechst_pos_G1S_rep1, log2fc_hoechst_pos_G1S_rep2, log2fc_hoechst_pos_G1S_mean = uti.get_average_delta_log2fc(data_group, data_cc, 'hoechst_pos_G1S_normalized')
        hoechst_pos_S_rep1, hoechst_pos_S_rep2, hoechst_pos_S_mean, _, _, _, log2fc_hoechst_pos_S_rep1, log2fc_hoechst_pos_S_rep2, log2fc_hoechst_pos_S_mean = uti.get_average_delta_log2fc(data_group, data_cc, 'hoechst_pos_S_normalized')
        hoechst_pos_G2M_rep1, hoechst_pos_G2M_rep2, hoechst_pos_G2M_mean, _, _, _, log2fc_hoechst_pos_G2M_rep1, log2fc_hoechst_pos_G2M_rep2, log2fc_hoechst_pos_G2M_mean = uti.get_average_delta_log2fc(data_group, data_cc, 'hoechst_pos_G2M_normalized')
        hoechst_pos_G2MG1_rep1, hoechst_pos_G2MG1_rep2, hoechst_pos_G2MG1_mean, _, _, _, log2fc_hoechst_pos_G2MG1_rep1, log2fc_hoechst_pos_G2MG1_rep2, log2fc_hoechst_pos_G2MG1_mean = uti.get_average_delta_log2fc(data_group, data_cc, 'hoechst_pos_G2MG1_normalized')

        df.loc[len(df.index)] = [batch, i+1, data_group['treatment'][0],
                                 hoechst_rep1, hoechst_rep2, hoechst_mean,
                                 log2fc_hoechst_rep1, log2fc_hoechst_rep2, log2fc_hoechst_mean,
                                 n_rep1, n_rep2, n_mean,
                                 log2fc_n_rep1, log2fc_n_rep2, log2fc_n_mean,
                                 n_neg_rep1, n_neg_rep2, n_neg_mean,
                                 log2fc_n_neg_rep1, log2fc_n_neg_rep2, log2fc_n_neg_mean,
                                 n_pos_rep1, n_pos_rep2, n_pos_mean,
                                 log2fc_n_pos_rep1, log2fc_n_pos_rep2, log2fc_n_pos_mean,
                                 ratio_rep1, ratio_rep2, ratio_mean,
                                 log2fc_ratio_rep1, log2fc_ratio_rep2, log2fc_ratio_mean,
                                 per_G1_rep1, per_G1_rep2, per_G1_mean,
                                 per_G1S_rep1, per_G1S_rep2, per_G1S_mean,
                                 per_S_rep1, per_S_rep2, per_S_mean,
                                 per_G2M_rep1, per_G2M_rep2, per_G2M_mean,
                                 per_G2MG1_rep1, per_G2MG1_rep2, per_G2MG1_mean,
                                 delta_G1_rep1, delta_G1_rep2, delta_G1_mean,
                                 delta_G1S_rep1, delta_G1S_rep2, delta_G1S_mean,
                                 delta_S_rep1, delta_S_rep2, delta_S_mean,
                                 delta_G2M_rep1, delta_G2M_rep2, delta_G2M_mean,
                                 delta_G2MG1_rep1, delta_G2MG1_rep2, delta_G2MG1_mean,
                                 per_neg_G1_rep1, per_neg_G1_rep2, per_neg_G1_mean,
                                 per_neg_G1S_rep1, per_neg_G1S_rep2, per_neg_G1S_mean,
                                 per_neg_S_rep1, per_neg_S_rep2, per_neg_S_mean,
                                 per_neg_G2M_rep1, per_neg_G2M_rep2, per_neg_G2M_mean,
                                 per_neg_G2MG1_rep1, per_neg_G2MG1_rep2, per_neg_G2MG1_mean,
                                 delta_neg_G1_rep1, delta_neg_G1_rep2, delta_neg_G1_mean,
                                 delta_neg_G1S_rep1, delta_neg_G1S_rep2, delta_neg_G1S_mean,
                                 delta_neg_S_rep1, delta_neg_S_rep2, delta_neg_S_mean,
                                 delta_neg_G2M_rep1, delta_neg_G2M_rep2, delta_neg_G2M_mean,
                                 delta_neg_G2MG1_rep1, delta_neg_G2MG1_rep2, delta_neg_G2MG1_mean,
                                 per_pos_G1_rep1, per_pos_G1_rep2, per_pos_G1_mean,
                                 per_pos_G1S_rep1, per_pos_G1S_rep2, per_pos_G1S_mean,
                                 per_pos_S_rep1, per_pos_S_rep2, per_pos_S_mean,
                                 per_pos_G2M_rep1, per_pos_G2M_rep2, per_pos_G2M_mean,
                                 per_pos_G2MG1_rep1, per_pos_G2MG1_rep2, per_pos_G2MG1_mean,
                                 delta_pos_G1_rep1, delta_pos_G1_rep2, delta_pos_G1_mean,
                                 delta_pos_G1S_rep1, delta_pos_G1S_rep2, delta_pos_G1S_mean,
                                 delta_pos_S_rep1, delta_pos_S_rep2, delta_pos_S_mean,
                                 delta_pos_G2M_rep1, delta_pos_G2M_rep2, delta_pos_G2M_mean,
                                 delta_pos_G2MG1_rep1, delta_pos_G2MG1_rep2, delta_pos_G2MG1_mean,
                                 delta_sum_rep1, delta_sum_rep2, delta_sum_mean,
                                 delta_neg_sum_rep1, delta_neg_sum_rep2, delta_neg_sum_mean,
                                 delta_pos_sum_rep1, delta_pos_sum_rep2, delta_pos_sum_mean,
                                 cc_filter,
                                 hoechst_G1_rep1, hoechst_G1_rep2, hoechst_G1_mean,
                                 hoechst_G1S_rep1, hoechst_G1S_rep2, hoechst_G1S_mean,
                                 hoechst_S_rep1, hoechst_S_rep2, hoechst_S_mean,
                                 hoechst_G2M_rep1, hoechst_G2M_rep2, hoechst_G2M_mean,
                                 hoechst_G2MG1_rep1, hoechst_G2MG1_rep2, hoechst_G2MG1_mean,
                                 hoechst_neg_G1_rep1, hoechst_neg_G1_rep2, hoechst_neg_G1_mean,
                                 hoechst_neg_G1S_rep1, hoechst_neg_G1S_rep2, hoechst_neg_G1S_mean,
                                 hoechst_neg_S_rep1, hoechst_neg_S_rep2, hoechst_neg_S_mean,
                                 hoechst_neg_G2M_rep1, hoechst_neg_G2M_rep2, hoechst_neg_G2M_mean,
                                 hoechst_neg_G2MG1_rep1, hoechst_neg_G2MG1_rep2, hoechst_neg_G2MG1_mean,
                                 hoechst_pos_G1_rep1, hoechst_pos_G1_rep2, hoechst_pos_G1_mean,
                                 hoechst_pos_G1S_rep1, hoechst_pos_G1S_rep2, hoechst_pos_G1S_mean,
                                 hoechst_pos_S_rep1, hoechst_pos_S_rep2, hoechst_pos_S_mean,
                                 hoechst_pos_G2M_rep1, hoechst_pos_G2M_rep2, hoechst_pos_G2M_mean,
                                 hoechst_pos_G2MG1_rep1, hoechst_pos_G2MG1_rep2, hoechst_pos_G2MG1_mean,
                                 log2fc_hoechst_G1_rep1, log2fc_hoechst_G1_rep2, log2fc_hoechst_G1_mean,
                                 log2fc_hoechst_G1S_rep1, log2fc_hoechst_G1S_rep2, log2fc_hoechst_G1S_mean,
                                 log2fc_hoechst_S_rep1, log2fc_hoechst_S_rep2, log2fc_hoechst_S_mean,
                                 log2fc_hoechst_G2M_rep1, log2fc_hoechst_G2M_rep2, log2fc_hoechst_G2M_mean,
                                 log2fc_hoechst_G2MG1_rep1, log2fc_hoechst_G2MG1_rep2, log2fc_hoechst_G2MG1_mean,
                                 log2fc_hoechst_neg_G1_rep1, log2fc_hoechst_neg_G1_rep2, log2fc_hoechst_neg_G1_mean,
                                 log2fc_hoechst_neg_G1S_rep1, log2fc_hoechst_neg_G1S_rep2, log2fc_hoechst_neg_G1S_mean,
                                 log2fc_hoechst_neg_S_rep1, log2fc_hoechst_neg_S_rep2, log2fc_hoechst_neg_S_mean,
                                 log2fc_hoechst_neg_G2M_rep1, log2fc_hoechst_neg_G2M_rep2, log2fc_hoechst_neg_G2M_mean,
                                 log2fc_hoechst_neg_G2MG1_rep1, log2fc_hoechst_neg_G2MG1_rep2, log2fc_hoechst_neg_G2MG1_mean,
                                 log2fc_hoechst_pos_G1_rep1, log2fc_hoechst_pos_G1_rep2, log2fc_hoechst_pos_G1_mean,
                                 log2fc_hoechst_pos_G1S_rep1, log2fc_hoechst_pos_G1S_rep2, log2fc_hoechst_pos_G1S_mean,
                                 log2fc_hoechst_pos_S_rep1, log2fc_hoechst_pos_S_rep2, log2fc_hoechst_pos_S_mean,
                                 log2fc_hoechst_pos_G2M_rep1, log2fc_hoechst_pos_G2M_rep2, log2fc_hoechst_pos_G2M_mean,
                                 log2fc_hoechst_pos_G2MG1_rep1, log2fc_hoechst_pos_G2MG1_rep2, log2fc_hoechst_pos_G2MG1_mean
                                 ]
    else:
        df.loc[len(df.index)] = [batch, i+1, data_group['treatment'][0]] + ['NA'] * 129 + [0] + ['NA'] * 90

df = df.copy()
df['qc_filter'] = [0 if df['hoechst_rep1'][i] == 'NA' else 1 for i in range(len(df))]
print(len(df))
print(len(df[df['cc_filter'] == 0]))
print(df[df['cc_filter'] == 0]['treatment'])
print(len(df[df['hoechst_rep1'] == 'NA']))

df.to_csv('%s/01_summary/%s_grouping.txt' % (output_dir, batch), index=False, sep='\t')
print("DONE!")