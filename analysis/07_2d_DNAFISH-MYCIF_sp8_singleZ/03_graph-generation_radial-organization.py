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
           'radial_curve_nuclear', 'radial_curve_DNAFISH', 'angle_curve_nuclear', 'angle_curve_DNAFISH',
           'percentage_area_curve_ecDNA', 'percentage_area_ratio_curve_ecDNA', 'percentage_int_curve_ecDNA',
           'cum_area_ind_ecDNA', 'cum_area_ratio_ind_ecDNA', 'cum_int_ind_ecDNA', 'cum_area_ind_ecDNA_filled',
           'cum_area_ratio_ind_ecDNA_filled', 'cum_int_ind_ecDNA_filled']

# LOAD FILE
data = pd.read_csv('%s%s_v%s.txt' % (master_folder, sample, version), na_values=['.'], sep='\t')

for f in feature:
    data[f] = [dat.str_to_float(data[f][i]) for i in range(len(data))]

# test
plt.subplots(figsize=(9, 6))
plt.scatter(data['total_area_ecDNA'], data['dis_to_hub_area'], c=data['radial_edge'])
plt.show()

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

print("DONE!")