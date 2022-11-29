import pandas as pd
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt
import seaborn as sns
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221017_periphery-localization_analysis/20220826_neuroblastoma/A/"
"""sample_lst = ['AU', 'AX', 'BT', 'BZ', 'C', 'CF', 'CN', 'CQ', 'DK', 'DM', 'DY', 'DZ', 'EF', 'EH', 'EM', 'ER', 'F', 'FI',
              'FL', 'FP', 'FV', 'FW', 'FX', 'GA', 'GY', 'GZ', 'HC', 'HL', 'HM', 'IP', 'IU', 'IX', 'JZ', 'KI', 'KJ',
              'KV', 'LM', 'ME', 'MU']"""
# sample_lst = ['AC', 'AK', 'CB', 'EB', 'FR', 'GL', 'GO', 'IC', 'NE', 'V']
group = 'A'
version = 1
save_folder = '%sfig/' % master_folder

# default setting
hist_colors = [(0.95, 0.50, 0.50), (0.90, 0.90, 0.90), (0.50, 0.90, 0.90)]
line_colors = [(0.85, 0.35, 0.25), (0.30, 0.30, 0.30), (0.30, 0.70, 0.70)]
line_color_4 = [(0.85, 0.35, 0.25), (0.30, 0.30, 0.30), (0.30, 0.70, 0.70), (0.95, 0.85, 0)]
# '#FFA500', '#40E0D0', '#DA70D6', '#00CED1'

feature = ['mean_int_ind_ecDNA', 'total_int_ind_ecDNA', 'area_ind_ecDNA', 'area_ratio_ind_ecDNA', 'g',
           'radial_curve_nuclear', 'radial_curve_DNAFISH', 'angle_curve_nuclear', 'angle_curve_DNAFISH',
           'percentage_area_curve_ecDNA', 'percentage_area_ratio_curve_ecDNA', 'percentage_int_curve_ecDNA',
           'cum_area_ind_ecDNA', 'cum_area_ratio_ind_ecDNA', 'cum_int_ind_ecDNA', 'cum_area_ind_ecDNA_filled',
           'cum_area_ratio_ind_ecDNA_filled', 'cum_int_ind_ecDNA_filled']

# LOAD FILE
data_ori = pd.read_csv('%s%s_v%s.txt' % (master_folder, group, version), na_values=['.'], sep='\t')
data = data_ori[data_ori['area_nuclear'] > 100].copy().reset_index(drop=True)

for f in feature:
    data[f] = [dat.str_to_float(data[f][i]) for i in range(len(data))]

data['radial_distance'] = data['radial_edge']-data['radial_center']
data['normalized_dis_to_hub_area_v2'] = data['dis_to_hub_area_v2']/data['normalized_r']
data['normalized_cluster'] = data['normalized_dis_to_hub_area_v2']/(data['total_area_ecDNA']**0.5)
data['radial_curve_normalized'] = \
    [list(np.array(data['radial_curve_DNAFISH'][i])/np.array(data['radial_curve_nuclear'][i])) for i in range(len(data))]

# radial curve
print("Plotting radial curve...")
x = np.arange(0.01, 0.99, 0.01)
x_label = 'relative r'
up_level = 90

number_nuclear = len(data)

mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data['radial_curve_DNAFISH'].tolist())
mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(data['radial_curve_nuclear'].tolist())
mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_normalized'].tolist())

plt.subplots(figsize=(12, 9))
for i in range(len(data)):
    plt.plot(x[:up_level], data['radial_curve_DNAFISH'][i][:up_level], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
for i in range(len(data)):
    plt.plot(x[:up_level], data['radial_curve_nuclear'][i][:up_level], alpha=0.05, color=[line_colors[1][j]+0.05 for j in range(len(line_colors[1]))])
plt.plot(x[:up_level], mean_curve1[:up_level], color=line_colors[0], label='%s, n=%s' % ('ecDNA (MYC DNA FISH)', number_nuclear))
plt.plot(x[:up_level], mean_curve2[:up_level], color=line_colors[1], label='%s, n=%s' % ('DNA (hoechst stain)', number_nuclear))
plt.plot(x[:up_level], ci_lower1[:up_level], color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x[:up_level], ci_higher1[:up_level], color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x[:up_level], ci_lower2[:up_level], color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x[:up_level], ci_higher2[:up_level], color=line_colors[1], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylim([0.75, 1.25])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s/%s_radial_curve.pdf' % (save_folder, group))
plt.show()

plt.subplots(figsize=(12, 9))
for i in range(len(data)):
    plt.plot(x[:up_level], data['radial_curve_normalized'][i][:up_level], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
plt.plot(x[:up_level], mean_curve3[:up_level], color=line_colors[0], label='%s, n=%s' % ('ecDNA (normalized)', number_nuclear))
plt.plot(x[:up_level], ci_lower3[:up_level], color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x[:up_level], ci_higher3[:up_level], color=line_colors[0], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.75, 1.25])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s/%s_radial_curve_normalized.pdf' % (save_folder, group))
plt.show()

# radial distribution vs cluster degree
"""plt.subplots(figsize=(12, 9))
plt.scatter(data['max_area_ecDNA'], data['radial_distance'], alpha=0.8)
plt.xlabel('max_area_ecDNA')
plt.ylabel('radial_distance')
plt.xlim([-5, 505])
plt.savefig('%s/%s_radial_distance-vs-max_area_ecDNA.pdf' % (save_folder, group))
plt.show()"""

data_filtered = data[data['n_ecDNA'] > 1].copy().reset_index(drop=True)

plt.subplots(figsize=(12, 9))
plt.scatter(data_filtered['total_area_ratio_ecDNA']**0.5, data_filtered['normalized_dis_to_hub_area_v2'],
            c=data_filtered['max_area_ecDNA'], vmin=0, vmax=500)
plt.xlabel('total_area_ratio_ecDNA_sqrt')
plt.ylabel('normalized_dis_to_hub_area_v2')
plt.ylim([0, 4])
plt.xlim([0, 1])
plt.savefig('%s/%s_dis_to_hub_area_v2-vs-total_area_ratio_ecDNA_sqrt-color-max_area_ecDNA.pdf' % (save_folder, group))
plt.show()

plt.subplots(figsize=(12, 9))
plt.scatter(data_filtered['radial_distance'], data_filtered['normalized_cluster'],
            c=data_filtered['max_area_ecDNA'], vmin=0, vmax=500)
plt.xlabel('radial_distance')
plt.ylabel('normalized_cluster')
plt.savefig('%s/%s_radial_distance-vs-normalized_cluster.pdf' % (save_folder, group))
plt.show()

print("DONE!")