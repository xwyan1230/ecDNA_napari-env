import pandas as pd
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from pandas.plotting import parallel_coordinates
import shared.display as dis
import random
import shared.math as mat
import math
import seaborn as sns
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220916_double_funnel_test/"
sample = 'slide1_B2'
sample_prefix = ''
version = 1
save_folder = "%sv1_figures/%s/" % (master_folder, sample)
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

# default setting
hist_colors = [(0.95, 0.50, 0.50), (0.90, 0.90, 0.90), (0.50, 0.90, 0.90)]
line_colors = [(0.85, 0.35, 0.25), (0.30, 0.30, 0.30), (0.30, 0.70, 0.70)]
line_color_4 = [(0.85, 0.35, 0.25), (0.30, 0.30, 0.30), (0.30, 0.70, 0.70), (0.95, 0.85, 0)]
# '#FFA500', '#40E0D0', '#DA70D6', '#00CED1'
feature = ['mean_int_ind_ecDNA', 'total_int_ind_ecDNA', 'area_ind_ecDNA', 'area_ratio_ind_ecDNA', 'g',
           'radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_nuclear_20', 'radial_curve_DNAFISH_20',
           'angle_curve_nuclear', 'angle_curve_DNAFISH',
           'percentage_area_curve_ecDNA', 'percentage_area_ratio_curve_ecDNA', 'percentage_int_curve_ecDNA',
           'cum_area_ind_ecDNA', 'cum_area_ratio_ind_ecDNA', 'cum_int_ind_ecDNA', 'cum_area_ind_ecDNA_filled',
           'cum_area_ratio_ind_ecDNA_filled', 'cum_int_ind_ecDNA_filled']

# LOAD FILE
data = pd.read_csv('%s%s_v%s.txt' % (master_folder, sample, version), na_values=['.'], sep='\t')

for f in feature:
    data[f] = [dat.str_to_float(data[f][i]) for i in range(len(data))]

data['radial_distance'] = data['radial_edge']-data['radial_center']
data['normalized_r'] = (data['area_nuclear']/math.pi)**0.5

plt.subplots(figsize=(9, 6))
plt.scatter(data['dis_to_hub_area'], data['dis_to_hub_area_v2'])
plt.xlabel('dis_to_hub_area')
plt.ylabel('dis_to_hub_area_v2')
plt.savefig('%s/dis_to_hub_area_algorithm-update_%s.pdf' % (save_folder, sample))
plt.close()

# radial distribution vs cluster degree

plt.subplots(figsize=(9, 6))
plt.scatter(data['total_area_ecDNA'], data['dis_to_hub_area_v2'], c=data['max_area_ecDNA'], vmin=0, vmax=2500)
plt.xlabel('total_area_ecDNA')
plt.ylabel('dis_to_hub_area_v2')
# plt.savefig('%s/dis_to_hub_area_v2-vs-total_area_ecDNA-color-max_area_ecDNA_%s.pdf' % (save_folder, sample))
plt.close()

plt.subplots(figsize=(9, 6))
plt.scatter(data['dis_to_hub_area_v2']/data['total_area_ecDNA'], data['radial_distance'], c=data['max_area_ecDNA'], vmin=0, vmax=2500)
plt.xlabel('dis_to_hub_area_v2/total_area_ecDNA')
plt.ylabel('radial_distance')
plt.xlim([0, 0.15])
# plt.savefig('%s/radial_distance-vs-dis_to_hub_area_v2_normalized_by_total_area_ecDNA-color-max_area_ecDNA_%s.pdf' % (save_folder, sample))
plt.close()

plt.subplots(figsize=(9, 6))
plt.scatter((data['dis_to_hub_area_v2']/data['normalized_r'])/(data['total_area_ratio_ecDNA']**0.5), data['radial_distance'], c=data['max_area_ecDNA'], vmin=0, vmax=2500)
plt.xlabel('normalized_cluster')
plt.ylabel('radial_distance')
plt.savefig('%s/radial_distance-vs-normalized_cluster-color-max_area_ecDNA_%s.pdf' % (save_folder, sample))
plt.close()

plt.subplots(figsize=(9, 6))
plt.scatter(data['max_area_ecDNA'], data['radial_distance'], alpha=0.7)
plt.xlabel('max_area_ecDNA')
plt.ylabel('radial_distance')
# plt.savefig('%s/radial_distance-vs-max_area_ecDNA_%s.pdf' % (save_folder, sample))
plt.close()

# radial curve
print("Plotting radial curve...")

x = np.arange(0.01, 0.99, 0.01)
x_label = 'relative r'

number_nuclear = len(data)
data['radial_curve_normalized'] = \
    [list(np.array(data['radial_curve_DNAFISH'][i])/np.array(data['radial_curve_nuclear'][i])) for i in range(len(data))]

mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data['radial_curve_DNAFISH'].tolist())
mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(data['radial_curve_nuclear'].tolist())
mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_normalized'].tolist())

plt.subplots(figsize=(9, 6))
for i in range(len(data)):
    plt.plot(x, data['radial_curve_DNAFISH'][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
for i in range(len(data)):
    plt.plot(x, data['radial_curve_nuclear'][i], alpha=0.05, color=[line_colors[1][j]+0.05 for j in range(len(line_colors[1]))])
plt.plot(x, mean_curve1, color=line_colors[0], label='%s, n=%s' % ('ecDNA (MYC DNA FISH)', number_nuclear))
plt.plot(x, mean_curve2, color=line_colors[1], label='%s, n=%s' % ('DNA (hoechst stain)', number_nuclear))
plt.plot(x, ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x, ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylim([0, 2])
plt.ylabel('radial_curve')
plt.legend()
# plt.savefig('%s/radial_curve_%s.pdf' % (save_folder, sample))
plt.close()

plt.subplots(figsize=(9, 6))
for i in range(len(data)):
    plt.plot(x, data['radial_curve_normalized'][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
plt.plot(x, mean_curve3, color=line_colors[0], label='%s, n=%s' % ('ecDNA (normalized)', number_nuclear))
plt.plot(x, ci_lower3, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher3, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0, 2])
plt.ylabel('radial_curve')
plt.legend()
# plt.savefig('%s/radial_curve_normalized_%s.pdf' % (save_folder, sample))
plt.close()

x = np.arange(0, 1, 0.05)
x_label = 'relative r'

number_nuclear = len(data)
data['radial_curve_normalized_20'] = \
    [list(np.array(data['radial_curve_DNAFISH_20'][i])/np.array(data['radial_curve_nuclear_20'][i])) for i in range(len(data))]

mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data['radial_curve_DNAFISH_20'].tolist())
mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(data['radial_curve_nuclear_20'].tolist())
mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_normalized_20'].tolist())

plt.subplots(figsize=(9, 6))
for i in range(len(data)):
    plt.plot(x, data['radial_curve_DNAFISH_20'][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
for i in range(len(data)):
    plt.plot(x, data['radial_curve_nuclear_20'][i], alpha=0.05, color=[line_colors[1][j]+0.05 for j in range(len(line_colors[1]))])
plt.plot(x, mean_curve1, color=line_colors[0], label='%s, n=%s' % ('ecDNA (MYC DNA FISH)', number_nuclear))
plt.plot(x, mean_curve2, color=line_colors[1], label='%s, n=%s' % ('DNA (hoechst stain)', number_nuclear))
plt.plot(x, ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x, ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylim([0, 2])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s/radial_curve_20_%s.pdf' % (save_folder, sample))
plt.close()

plt.subplots(figsize=(9, 6))
for i in range(len(data)):
    plt.plot(x, data['radial_curve_normalized_20'][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
plt.plot(x, mean_curve3, color=line_colors[0], label='%s, n=%s' % ('ecDNA (normalized)', number_nuclear))
plt.plot(x, ci_lower3, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher3, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0, 2])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s/radial_curve_normalized_20_%s.pdf' % (save_folder, sample))
plt.close()

data_heatmap = pd.DataFrame(columns=['0-5', '5-10', '10-15', '15-20', '20-25', '25-30', '30-35', '35-40', '40-45',
                                     '45-50', '50-55', '55-60', '60-65', '65-70', '70-75', '75-80', '80-85',
                                     '85-90', '90-95', '95-100'])
data_heatmap.loc[0] = mean_curve3

plt.subplots(figsize=(12, 9))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=1.2, vmin=0.8, square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_%s.pdf' % (save_folder, sample))
plt.close()

print("DONE!")