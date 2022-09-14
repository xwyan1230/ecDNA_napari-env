import numpy as np
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
import skimage.io as skio
import shared.display as dis
import shared.dataframe as dat
import shared.image as ima
import pandas as pd
import napari
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220826_neuroblastoma/"
group = 'C'
save_path = master_folder
version = 1

data = pd.read_csv('%s%s_v%s.txt' % (master_folder, group, version), na_values=['.'], sep='\t')

sample_lst = list(set(data['sample'].tolist()))
sample_size = len(sample_lst)

feature_str = ['radial_curve_DNAFISH', 'radial_curve_nuclear', 'cum_area_ind_ecDNA', 'cum_area_ratio_ind_ecDNA', 'cum_area_ind_ecDNA_filled',
               'cum_area_ratio_ind_ecDNA_filled', 'percentage_area_curve_ecDNA', 'percentage_area_ratio_curve_ecDNA',
               'area_ind_ecDNA', 'area_ratio_ind_ecDNA', 'angle_curve_nuclear', 'angle_curve_DNAFISH']
for f in feature_str:
    data[f] = [dat.str_to_float(data[f][i]) for i in range(len(data))]
data['radial_curve'] = [list(np.array(data['radial_curve_DNAFISH'].tolist()[i])/np.array(data['radial_curve_nuclear'].tolist()[i])) for i in range(len(data))]

data_mean = pd.DataFrame(columns=['sample', 'n_nuclear', 'area_nuclear', 'n_ecDNA',
                                  'total_area_ecDNA', 'total_area_ratio_ecDNA', 'area_ind_ecDNA', 'area_ratio_ind_ecDNA',
                                  'max_area_ecDNA', 'max_area_ratio_ecDNA',
                                  'radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve',
                                  'radial_center', 'radial_edge', 'relative_r_area',
                                  'angle_curve_nuclear', 'angle_curve_DNAFISH', 'angle_value',
                                  'dis_to_hub_area',
                                  'percentage_area_curve_ecDNA', 'percentage_area_n_half',
                                  'percentage_area_ratio_curve_ecDNA', 'percentage_area_ratio_n_half',
                                  'cum_area_ind_ecDNA', 'cum_area_n_half',
                                  'cum_area_ratio_ind_ecDNA', 'cum_area_ratio_n_half'])
for sample in sample_lst:
    data_sample = data[data['sample'] == sample].copy().reset_index()
    if len(data_sample) != 0:
        n_nuclear = len(data_sample)
        area_nuclear = np.sum(data_sample['area_nuclear'].tolist())/n_nuclear
        n_ecDNA = np.sum(data_sample['n_ecDNA'].tolist())/n_nuclear
        total_area_ecDNA = np.sum(data_sample['total_area_ecDNA'].tolist())/n_nuclear
        total_area_ratio_ecDNA = np.sum(data_sample['total_area_ratio_ecDNA'].tolist()) / n_nuclear
        max_area_ecDNA = np.sum(data_sample['max_area_ecDNA'].tolist())/n_nuclear
        max_area_ratio_ecDNA = np.sum(data_sample['max_area_ratio_ecDNA'].tolist())/n_nuclear
        radial_center = np.sum(data_sample['radial_center'].tolist()) / n_nuclear
        radial_edge = np.sum(data_sample['radial_edge'].tolist()) / n_nuclear
        relative_r_area = np.sum(data_sample['relative_r_area'].tolist()) / n_nuclear
        angle_value = np.sum(data_sample['angle_value'].tolist()) / n_nuclear
        dis_to_hub_area = np.sum(data_sample['dis_to_hub_area'].tolist()) / n_nuclear
        percentage_area_n_half = np.sum(data_sample['percentage_area_n_half'].tolist()) / n_nuclear
        percentage_area_ratio_n_half = np.sum(data_sample['percentage_area_ratio_n_half'].tolist()) / n_nuclear
        cum_area_n_half = np.sum(data_sample['cum_area_n_half'].tolist()) / n_nuclear
        cum_area_ratio_n_half = np.sum(data_sample['cum_area_ratio_n_half'].tolist()) / n_nuclear

        radial_curve_nuclear, _, _ = dat.mean_list(data_sample['radial_curve_nuclear'].tolist())
        radial_curve_DNAFISH, _, _ = dat.mean_list(data_sample['radial_curve_DNAFISH'].tolist())
        radial_curve, _, _ = dat.mean_list(data_sample['radial_curve'].tolist())
        angle_curve_nuclear, _, _ = dat.mean_list(data_sample['angle_curve_nuclear'].tolist())
        angle_curve_DNAFISH, _, _ = dat.mean_list(data_sample['angle_curve_DNAFISH'].tolist())
        percentage_area_curve_ecDNA, _, _ = dat.mean_list(dat.list_fill_with_last_num(data_sample['percentage_area_curve_ecDNA'].tolist()))
        percentage_area_ratio_curve_ecDNA, _, _ = dat.mean_list(dat.list_fill_with_last_num(data_sample['percentage_area_ratio_curve_ecDNA'].tolist()))
        cum_area_ind_ecDNA, _, _ = dat.mean_list(dat.list_fill_with_last_num(data_sample['cum_area_ind_ecDNA'].tolist()))
        cum_area_ratio_ind_ecDNA, _, _ = dat.mean_list(dat.list_fill_with_last_num(data_sample['cum_area_ratio_ind_ecDNA'].tolist()))

        area_ind_ecDNA = []
        area_ind_ecDNA = [area_ind_ecDNA + data_sample['area_ind_ecDNA'].tolist()[i] for i in range(len(data_sample))]
        area_ratio_ind_ecDNA = []
        area_ratio_ind_ecDNA = [area_ratio_ind_ecDNA + data_sample['area_ratio_ind_ecDNA'].tolist()[i] for i in range(len(data_sample))]

        data_mean.loc[len(data_mean.index)] = [sample, n_nuclear, area_nuclear, n_ecDNA,
                                               total_area_ecDNA, total_area_ratio_ecDNA, area_ind_ecDNA, area_ratio_ind_ecDNA,
                                               max_area_ecDNA, max_area_ratio_ecDNA,
                                               radial_curve_nuclear, radial_curve_DNAFISH, radial_curve,
                                               radial_center, radial_edge, relative_r_area,
                                               angle_curve_nuclear, angle_curve_DNAFISH, angle_value,
                                               dis_to_hub_area, percentage_area_curve_ecDNA, percentage_area_n_half,
                                               percentage_area_ratio_curve_ecDNA, percentage_area_ratio_n_half,
                                               cum_area_ind_ecDNA, cum_area_n_half,
                                               cum_area_ratio_ind_ecDNA, cum_area_ratio_n_half]

data_mean['group'] = [group] * len(data_mean)
data_mean.to_csv('%s%s_v%s_mean.txt' % (master_folder, group, version), index=False, sep='\t')

print("DONE!")