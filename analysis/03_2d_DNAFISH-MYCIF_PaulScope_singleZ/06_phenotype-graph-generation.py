import pandas as pd
import shared.display as dis
import numpy as np
import os

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220726_BRDfamily_screen/"
features = ['area_nuclear', 'mean_int_nuclear', 'total_int_nuclear', 'mean_int_IF', 'total_int_IF', 'n_ecDNA',
            'mean_int_DNAFISH_norm', 'total_int_DNAFISH_norm', 'mean_int_ecDNA_norm', 'total_int_ecDNA_norm',
            'area_ratio_ecDNA', 'max_area_ratio_ecDNA', 'radial_center', 'radial_edge', 'percentage_area_n_half',
            'percentage_area_ratio_n_half', 'percentage_int_n_half', 'percentage_int_norm_n_half', 'cum_area_n_half',
            'cum_area_ratio_n_half', 'cum_int_n_half', 'cum_int_norm_n_half', 'dis_to_hub_area', 'dis_to_hub_int',
            'dis_to_hub_int_norm', 'mean_int_ind_ecDNA_norm', 'total_int_ind_ecDNA_norm', 'area_ratio_ind_ecDNA',
            'g_value', 'angle_value', 'relative_r_area', 'relative_r_int', 'relative_r_int_norm']
save_folder = "%sv2_volcano/" % master_folder
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

data_p_full = pd.read_csv(("%sv2_summary/summary_p.txt" % master_folder), na_values=['.'], sep='\t')
data_mean_full = pd.read_csv(("%sv2_summary/summary_mean.txt" % master_folder), na_values=['.'], sep='\t')
data_gamma_full = pd.read_csv(("%sv2_summary/summary_gamma.txt" % master_folder), na_values=['.'], sep='\t')

data_p = data_p_full[(data_p_full['gene'] != 'NTC1') & (data_p_full['gene'] != 'NTC2')]
data_mean = data_mean_full[(data_mean_full['gene'] != 'NTC1') & (data_mean_full['gene'] != 'NTC2')]
data_gamma = data_gamma_full[(data_gamma_full['gene'] != 'NTC1') & (data_gamma_full['gene'] != 'NTC2')]

centers = []
for i in range(len(features)):
    centers.append(data_mean[features[i]].tolist()[-1])

data_value = data_mean[:-1].copy()

for i in range(len(features)):
    dis.plot_volcano_hit(data_p, data_value, data_gamma, features[i], centers[i], 30, 'Y', save_folder)

print("DONE!")
