import pandas as pd
import numpy as np
import shared.dataframe as dat
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from pandas.plotting import parallel_coordinates
import shared.display as dis
import random
import shared.math as mat
import seaborn as sns
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220726_BRDfamily_screen/"
version = 1
sample_row_lst = ['B', 'C', 'D', 'E', 'F', 'G']
sample_col_lst = np.arange(2, 12, 1)
sample_lst = []
for sample_row in sample_row_lst:
    for sample_col in sample_col_lst:
        sample_lst.append("%s%s" % (sample_row, sample_col))
save_folder = "%sv2_summary/" % master_folder
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

# default setting
control = ['C7', 'D6', 'E7', 'F6', 'G3', 'G6', 'G7', 'G10']
WT = ['B2', 'F7', 'E6', 'G11', 'B11', 'C6', 'D7', 'G2']
data_gene = pd.read_csv('%sgene.txt' % master_folder, na_values=['.'], sep='\t')
feature = ['radial_curve_DNAFISH', 'radial_curve_nuclear', 'angle_curve_DNAFISH', 'angle_curve_nuclear',
           'area_individual_ecDNA', 'area_ratio_individual_ecDNA', 'mean_int_individual_ecDNA',
           'mean_int_individual_ecDNA_norm', 'total_int_individual_ecDNA', 'total_int_individual_ecDNA_norm',
           'percentage_area_curve_ecDNA', 'percentage_area_ratio_curve_ecDNA', 'percentage_int_curve_ecDNA',
           'percentage_int_curve_ecDNA_norm', 'cum_area_ind_ecDNA', 'cum_area_ind_ecDNA_filled',
           'cum_area_ratio_ind_ecDNA', 'cum_area_ratio_ind_ecDNA_filled', 'cum_int_ind_ecDNA',
           'cum_int_ind_ecDNA_filled', 'cum_int_ind_ecDNA_norm', 'cum_int_ind_ecDNA_norm_filled',
           'angle_curve_DNAFISH_bg_correct', 'radial_curve_DNAFISH_bg_correct', 'g']
single_feature_lst = ['z_ratio', 'limit', 'bg_int', 'area_nuclear', 'mean_int_nuclear', 'total_int_nuclear',
                      'mean_int_IF', 'total_int_IF', 'n_ecDNA', 'mean_int_DNAFISH', 'mean_int_DNAFISH_norm',
                      'total_int_DNAFISH', 'total_int_DNAFISH_norm', 'mean_int_ecDNA', 'mean_int_ecDNA_norm',
                      'total_int_ecDNA', 'total_int_ecDNA_norm',  'total_area_ecDNA', 'area_ratio_ecDNA',
                      'max_area_ecDNA', 'max_area_ratio_ecDNA', 'radial_center', 'radial_edge',
                      'percentage_area_n_half', 'percentage_area_ratio_n_half', 'percentage_int_n_half',
                      'percentage_int_norm_n_half', 'cum_area_n_half', 'cum_area_ratio_n_half', 'cum_int_n_half',
                      'cum_int_norm_n_half', 'dis_to_hub_area', 'dis_to_hub_int', 'dis_to_hub_int_norm',
                      'g_value', 'angle_value', 'relative_r_area', 'relative_r_int', 'relative_r_int_norm']

multi_feature_lst = ['mean_int_ind_ecDNA', 'mean_int_ind_ecDNA_norm', 'total_int_ind_ecDNA', 'total_int_ind_ecDNA_norm',
                     'area_ind_ecDNA', 'area_ratio_ind_ecDNA']

data = pd.DataFrame()

data_WT = pd.DataFrame()
for i in range(len(WT)):
    data_temp = pd.read_csv('%s%s/%s/%s_v%s.txt' % (master_folder, WT[i][0], WT[i][1:], WT[i], version), na_values=['.'],
                            sep='\t')
    if len(data_temp) >= 100:
        data_WT = pd.concat([data_WT, data_temp.sample(n=100)], axis=0, ignore_index=True)
    else:
        data_WT = pd.concat([data_WT, data_temp], axis=0, ignore_index=True)
data_WT['sample'] = ['WT'] * len(data_WT)
data_WT['gene'] = ['WT'] * len(data_WT)
data = pd.concat([data, data_WT], axis=0, ignore_index=True)
for f in feature:
    data_WT[f] = [dat.str_to_float(data_WT[f][i]) for i in range(len(data_WT))]

data_WT_area_ind_ecDNA = dat.list_addup_from_df(data_WT, 'area_individual_ecDNA')
data_WT_mean_int_ind_ecDNA = dat.list_addup_from_df(data_WT, 'mean_int_individual_ecDNA')
data_WT_area_ratio_ind_ecDNA = dat.list_addup_from_df(data_WT, 'area_ratio_individual_ecDNA')
data_WT_mean_int_ind_ecDNA_norm = dat.list_addup_from_df(data_WT, 'mean_int_individual_ecDNA_norm')
data_WT_total_int_ind_ecDNA = dat.list_addup_from_df(data_WT, 'total_int_individual_ecDNA')
data_WT_total_int_ind_ecDNA_norm = dat.list_addup_from_df(data_WT, 'total_int_individual_ecDNA_norm')

data_m_WT = pd.DataFrame({'area_ind_ecDNA': data_WT_area_ind_ecDNA,
                          'area_ratio_ind_ecDNA': data_WT_area_ratio_ind_ecDNA,
                          'mean_int_ind_ecDNA': data_WT_mean_int_ind_ecDNA,
                          'mean_int_ind_ecDNA_norm': data_WT_mean_int_ind_ecDNA_norm,
                          'total_int_ind_ecDNA': data_WT_total_int_ind_ecDNA,
                          'total_int_ind_ecDNA_norm': data_WT_total_int_ind_ecDNA_norm,
                          'sample': ['WT'] * len(data_WT_area_ind_ecDNA)})

column_lst = ['sample', 'n_s', 'n_m'] + single_feature_lst + multi_feature_lst + ['gene']

data_mean = pd.DataFrame(columns=column_lst)
data_p_max = pd.DataFrame(columns=column_lst)
data_p = pd.DataFrame(columns=column_lst)
data_gamma = pd.DataFrame(columns=column_lst)

single_feature_WT_mean = []
for f in single_feature_lst:
    single_feature_WT_mean.append(np.mean(data_WT[f]))

multi_feature_WT_mean = []
for f in multi_feature_lst:
    multi_feature_WT_mean.append(np.mean(data_m_WT[f]))

for sample in sample_lst:
    print("Analyzing sample %s..." % sample)
    gene = data_gene[data_gene['sample'] == sample]['gene'].tolist()[0]
    data_sample = pd.read_csv('%s%s/%s/%s_v%s.txt' % (master_folder, sample[0], sample[1:], sample, version), na_values=['.'],
                              sep='\t')
    if len(data_sample) != 0:
        data_sample['sample'] = [sample] * len(data_sample)
        data_sample['gene'] = [gene] * len(data_sample)
        data = pd.concat([data, data_sample], axis=0, ignore_index=True)

        for f in feature:
            data_sample[f] = [dat.str_to_float(data_sample[f][i]) for i in range(len(data_sample))]

        data_sample_area_ind_ecDNA = dat.list_addup_from_df(data_sample, 'area_individual_ecDNA')
        data_sample_mean_int_ind_ecDNA = dat.list_addup_from_df(data_sample, 'mean_int_individual_ecDNA')
        data_sample_area_ratio_ind_ecDNA = dat.list_addup_from_df(data_sample, 'area_ratio_individual_ecDNA')
        data_sample_mean_int_ind_ecDNA_norm = dat.list_addup_from_df(data_sample, 'mean_int_individual_ecDNA_norm')
        data_sample_total_int_ind_ecDNA = dat.list_addup_from_df(data_sample, 'total_int_individual_ecDNA')
        data_sample_total_int_ind_ecDNA_norm = dat.list_addup_from_df(data_sample, 'total_int_individual_ecDNA_norm')

        data_m_sample = pd.DataFrame({'area_ind_ecDNA': data_sample_area_ind_ecDNA,
                                      'area_ratio_ind_ecDNA': data_sample_area_ratio_ind_ecDNA,
                                      'mean_int_ind_ecDNA': data_sample_mean_int_ind_ecDNA,
                                      'mean_int_ind_ecDNA_norm': data_sample_mean_int_ind_ecDNA_norm,
                                      'total_int_ind_ecDNA': data_sample_total_int_ind_ecDNA,
                                      'total_int_ind_ecDNA_norm': data_sample_total_int_ind_ecDNA_norm,
                                      'sample': [sample] * len(data_sample_area_ind_ecDNA)})

        single_feature_mean = []
        single_feature_p_100 = []
        single_feature_max_p = []
        single_feature_gamma = []
        for f in single_feature_lst:
            single_feature_mean.append(np.mean(data_sample[f]))
            p_max = -np.log(ks_2samp(data_sample[f].sample(n=len(data_sample)).tolist(),
                                     data_WT[f].sample(n=len(data_sample)).tolist())[1])
            single_feature_max_p.append(p_max)
            if len(data_sample) >= 100:
                p_100 = -np.log(ks_2samp(data_sample[f].sample(n=100).tolist(), data_WT[f].sample(n=100).tolist())[1])
                single_feature_p_100.append(p_100)
                single_feature_gamma.append(np.abs(np.mean(data_sample[f] - np.mean(data_WT[f]))) * p_100)
            else:
                single_feature_p_100.append(p_max)
                single_feature_gamma.append(np.abs(np.mean(data_sample[f] - np.mean(data_WT[f]))) * p_max)

        multi_feature_mean = []
        multi_feature_p_250 = []
        multi_feature_max_p = []
        multi_feature_gamma = []
        for f in multi_feature_lst:
            multi_feature_mean.append(np.mean(data_m_sample[f]))
            p_max = -np.log(ks_2samp(data_m_sample[f].sample(n=len(data_m_sample)).tolist(),
                                     data_m_WT[f].sample(n=len(data_m_sample)).tolist())[1])
            multi_feature_max_p.append(p_max)
            if len(data_m_sample) >= 250:
                p_250 = -np.log(ks_2samp(data_m_sample[f].sample(n=250).tolist(), data_m_WT[f].sample(n=250).tolist())[1])
                multi_feature_p_250.append(p_250)
                multi_feature_gamma.append(np.abs(np.mean(data_m_sample[f] - np.mean(data_m_WT[f]))) * p_250)
            else:
                multi_feature_p_250.append(p_max)
                multi_feature_gamma.append(np.abs(np.mean(data_m_sample[f] - np.mean(data_m_WT[f]))) * p_max)

        n_s = len(data_sample)
        n_m = len(data_m_sample)

        data_mean.loc[len(data_mean.index)] = [sample, n_s, n_m] + single_feature_mean + multi_feature_mean + [gene]
        data_p_max.loc[len(data_p_max.index)] = [sample, n_s, n_m] + single_feature_max_p + multi_feature_max_p + [gene]
        data_p.loc[len(data_p.index)] = [sample, n_s, n_m] + single_feature_p_100 + multi_feature_p_250 + [gene]
        data_gamma.loc[len(data_gamma.index)] = [sample, n_s, n_m] + single_feature_gamma + multi_feature_gamma + [gene]

data_mean.loc[len(data_mean.index)] = ['WT', len(data_WT), len(data_m_WT)] + single_feature_WT_mean + \
                                      multi_feature_WT_mean + ['WT']

data_mean.to_csv('%ssummary_mean.txt' % save_folder, index=False, sep='\t')
data_p_max.to_csv('%ssummary_p_max.txt' % save_folder, index=False, sep='\t')
data_p.to_csv('%ssummary_p.txt' % save_folder, index=False, sep='\t')
data_gamma.to_csv('%ssummary_gamma.txt' % save_folder, index=False, sep='\t')
data.to_csv('%ssummary.txt' % save_folder, index=False, sep='\t')

print("DONE!")
