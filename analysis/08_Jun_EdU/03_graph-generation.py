from skimage.measure import regionprops, label
import shared.image as ima
import skimage.io as skio
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt
import os
import math
import pandas as pd

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220927_Jun_EGFR_RPAs33p_Edu/EGFR_RPAs33p_Edu/"
sample1 = 'gbm39ec hu'
sample2 = 'gbm39ec con'
save_folder = "%sv1_figures/%s/" % (master_folder, 'gbm39')
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

# default setting
hist_colors = [(0.95, 0.50, 0.50), (0.90, 0.90, 0.90), (0.50, 0.90, 0.90)]
line_colors = [(0.85, 0.35, 0.25), (0.30, 0.30, 0.30), (0.30, 0.70, 0.70)]
line_color_4 = [(0.85, 0.35, 0.25), (0.30, 0.30, 0.30), (0.30, 0.70, 0.70), (0.95, 0.85, 0)]
# '#FFA500', '#40E0D0', '#DA70D6', '#00CED1'
feature = ['mean_int_ind_ecDNA', 'total_int_ind_ecDNA', 'area_ind_ecDNA', 'area_ratio_ind_ecDNA',
           'percentage_area_curve_ecDNA', 'percentage_area_ratio_curve_ecDNA', 'percentage_int_curve_ecDNA',
           'cum_area_ind_ecDNA', 'cum_area_ratio_ind_ecDNA', 'cum_int_ind_ecDNA', 'cum_area_ind_ecDNA_filled',
           'cum_area_ratio_ind_ecDNA_filled', 'cum_int_ind_ecDNA_filled']

# LOAD FILE
data1 = pd.read_csv('%s%s.txt' % (master_folder, sample1), na_values=['.'], sep='\t')
data2 = pd.read_csv('%s%s.txt' % (master_folder, sample2), na_values=['.'], sep='\t')

data1['normalized_r'] = (data1['area_nuclear']/math.pi)**0.5
data1['normalized_cluster'] = (data1['dis_to_hub_area_v2']/data1['normalized_r'])/(data1['total_area_ratio_ecDNA']**0.5)
data2['normalized_r'] = (data2['area_nuclear']/math.pi)**0.5
data2['normalized_cluster'] = (data2['dis_to_hub_area_v2']/data2['normalized_r'])/(data2['total_area_ratio_ecDNA']**0.5)

data1_G2 = data1[data1['mean_int_EdU'] < 8].copy()
data2_G2 = data2[data2['mean_int_EdU'] < 8].copy()

features = ['total_area_ecDNA', 'dis_to_hub_area', 'percentage_area_n_half', 'cum_area_n_half']

# plt.scatter(data1_G2['total_area_ecDNA'], data1_G2['dis_to_hub_area'], c=data1_G2['RPA_count'], cmap='Blues', s=8)
# plt.hist(data1['mean_int_EdU'], bins=50)
# plt.scatter(data2_G2['total_area_ecDNA'], data2_G2['dis_to_hub_area'], c=data2_G2['RPA_count'], cmap='Reds', s=8)
plt.figure(figsize=(12, 9))
plt.scatter(data1_G2['total_area_ecDNA'], data1_G2['mean_int_RPA'], cmap='Blues', s=8, label='gbm39ec hu')
plt.scatter(data2_G2['total_area_ecDNA'], data2_G2['mean_int_RPA'], cmap='Reds', s=8, label='gbm39ec con')
plt.xlabel('total_area_ecDNA')
plt.ylabel('RPA_mean_int')
# plt.ylim([0,20])
plt.savefig('%s/%s_total_area-RPA_mean_int.pdf' % (save_folder, 'ec'))
plt.show()

plt.figure(figsize=(12, 9))
plt.scatter(data1_G2['normalized_cluster'], data1_G2['mean_int_RPA'], cmap='Blues', s=8, label='gbm39ec hu')
plt.scatter(data2_G2['normalized_cluster'], data2_G2['mean_int_RPA'], cmap='Reds', s=8, label='gbm39ec con')
plt.xlabel('normalized_cluster')
plt.ylabel('RPA_mean_int')
# plt.ylim([0,20])
plt.savefig('%s/%s_cluster-RPA_mean_int.pdf' % (save_folder, 'ec'))
plt.show()
