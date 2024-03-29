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
# master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220816_Natasha_ColoDM_reimage/"
master_folder = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220825_POLR3D/20220825_POLR3Dtest/'
sample = 'KO'
version = 1
save_folder = "%sv1_fig/%s_WT-Control_1/" % (master_folder, sample)
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

# default setting
control = ['Control_2']
WT = ['Control_1']
hist_colors = [(0.95, 0.50, 0.50), (0.90, 0.90, 0.90), (0.50, 0.90, 0.90)]
line_colors = [(0.85, 0.35, 0.25), (0.30, 0.30, 0.30), (0.30, 0.70, 0.70)]
line_color_4 = [(0.85, 0.35, 0.25), (0.30, 0.30, 0.30), (0.30, 0.70, 0.70), (0.95, 0.85, 0)]
# '#FFA500', '#40E0D0', '#DA70D6', '#00CED1'
feature = ['radial_curve_DNAFISH', 'radial_curve_nuclear', 'angle_curve_DNAFISH', 'angle_curve_nuclear',
           'area_ind_ecDNA', 'area_ratio_ind_ecDNA', 'mean_int_ind_ecDNA', 'total_int_ind_ecDNA',
           'percentage_area_curve_ecDNA', 'percentage_area_ratio_curve_ecDNA', 'percentage_int_curve_ecDNA',
           'cum_area_ind_ecDNA', 'cum_area_ind_ecDNA_filled',
           'cum_area_ratio_ind_ecDNA', 'cum_area_ratio_ind_ecDNA_filled', 'cum_int_ind_ecDNA',
           'cum_int_ind_ecDNA_filled', 'g']

# LOAD FILE
data_sample = pd.read_csv('%s%s_v%s.txt' % (master_folder, sample, version), na_values=['.'], sep='\t')
data_control1 = pd.read_csv('%s%s_v%s.txt' % (master_folder, control[0], version),
                            na_values=['.'], sep='\t')
data_WT = pd.read_csv('%s%s_v%s.txt' % (master_folder, WT[0], version), na_values=['.'], sep='\t')

data_sample['sample'] = [sample] * len(data_sample)
data_control1['sample'] = [control[0]] * len(data_control1)
data_WT['sample'] = ['WT'] * len(data_WT)
data_sample['sample_heatmap'] = [0] * len(data_sample)
data_WT['sample_heatmap'] = [1] * len(data_WT)
data_control1['sample_heatmap'] = [2] * len(data_control1)

for f in feature:
    data_sample[f] = [dat.str_to_float(data_sample[f][i]) for i in range(len(data_sample))]
    data_control1[f] = [dat.str_to_float(data_control1[f][i]) for i in range(len(data_control1))]
    data_WT[f] = [dat.str_to_float(data_WT[f][i]) for i in range(len(data_WT))]

# calculate
data_sample_area_ind_ecDNA = dat.list_addup_from_df(data_sample, 'area_ind_ecDNA')
data_control1_area_ind_ecDNA = dat.list_addup_from_df(data_control1, 'area_ind_ecDNA')
data_WT_area_ind_ecDNA = dat.list_addup_from_df(data_WT, 'area_ind_ecDNA')

data_sample_mean_int_ind_ecDNA = dat.list_addup_from_df(data_sample, 'mean_int_ind_ecDNA')
data_control1_mean_int_ind_ecDNA = dat.list_addup_from_df(data_control1, 'mean_int_ind_ecDNA')
data_WT_mean_int_ind_ecDNA = dat.list_addup_from_df(data_WT, 'mean_int_ind_ecDNA')

data_sample_area_ratio_ind_ecDNA = dat.list_addup_from_df(data_sample, 'area_ratio_ind_ecDNA')
data_control1_area_ratio_ind_ecDNA = dat.list_addup_from_df(data_control1, 'area_ratio_ind_ecDNA')
data_WT_area_ratio_ind_ecDNA = dat.list_addup_from_df(data_WT, 'area_ratio_ind_ecDNA')

# normalize dataset
data = pd.concat([data_sample, data_WT, data_control1], axis=0, ignore_index=True)
data_feature = data.copy()
all_feature = data.columns
feature = ['mean_int_ecDNA', 'area_ratio_ecDNA', 'radial_center', 'relative_r_area', 'n_ecDNA',
           'cum_area_ratio_n_half', 'dis_to_hub_area']
drop_feature = [i for i in all_feature if i not in feature]
data_feature = data_feature.drop(drop_feature, axis=1)

scaler = StandardScaler()
data_feature_scale = pd.DataFrame(scaler.fit_transform(data_feature))
data_feature_scale.columns = data_feature.columns
data_norm = pd.concat([data_feature_scale, data['sample_heatmap']], axis=1)

# heat map
plt.subplots(figsize=(6, 50))
ax1 = sns.heatmap(data_norm, cbar=0, linewidths=2, vmax=3, vmin=-2, square=True, cmap='viridis')
plt.savefig('%s/heatmap.pdf' % save_folder)
plt.close()

# pca
# 1. rescale the data to a [0,1] range
# 2. standardize the data to have a zero mean and unit standard deviation (currently using this one)
# data_feature_scale = (data_feature - data_feature.min())/(data_feature.max() - data_feature.min())
pca = sklearnPCA(n_components=2)
transformed = pd.DataFrame(pca.fit_transform(data_feature_scale))
eigenvalues = pca.explained_variance_
print(eigenvalues)

plt.scatter(transformed[data['sample_heatmap'] == 0][0], transformed[data['sample_heatmap'] == 0][1], label=sample,
            color=line_colors[0], alpha=0.5, s=10)
plt.scatter(transformed[data['sample_heatmap'] == 1][0], transformed[data['sample_heatmap'] == 1][1], label='WT',
            color=line_colors[1], alpha=0.5, s=10)

plt.legend()
plt.savefig('%s/pca.pdf' % save_folder)
plt.close()

# lda
"""lda = LDA(n_components=2)
lda_transformed = pd.DataFrame(lda.fit_transform(data_feature_scale, data['sample_heatmap']))

# Plot all three series
plt.scatter(lda_transformed[data['sample_heatmap'] == 0][0], lda_transformed[data['sample_heatmap'] == 0][1],
            label=sample, c=line_colors[0], alpha=0.5, s=10)
plt.scatter(lda_transformed[data['sample_heatmap'] == 1][0], lda_transformed[data['sample_heatmap'] == 1][1],
            label='WT', c=line_colors[1], alpha=0.5, s=10)

plt.legend()
plt.savefig('%s/lda.pdf' % save_folder)
plt.close()"""

# parallel coordinates
print("Plotting pc...")
# Select features to include in the plot
plot_feat = data_feature.columns
# Concat classes with the normalized data
data_norm = pd.concat([data_feature_scale[plot_feat], data['sample']], axis=1)
# Perform parallel coordinate plot
plt.subplots(figsize=(20, 4))
parallel_coordinates(data_norm, 'sample', color=line_color_4, alpha=0.3)
# Display legend and show plot
plt.legend(loc=3)
plt.savefig('%s/pc.pdf' % save_folder)
plt.close()

# single value feature
print("Plotting single value feature...")
feature = ['z_ratio', 'limit', 'area_nuclear', 'mean_int_DNAFISH', 'mean_int_nuclear',
           'total_int_DNAFISH', 'g_value', 'angle_value',
           'total_int_nuclear', 'radial_center', 'radial_edge', 'total_area_ecDNA', 'total_area_ratio_ecDNA',
           'mean_int_ecDNA', 'total_int_ecDNA',
           'max_area_ecDNA', 'max_area_ratio_ecDNA', 'n_ecDNA', 'percentage_area_n_half',
           'percentage_area_ratio_n_half', 'percentage_int_n_half', 'cum_area_n_half',
           'cum_area_ratio_n_half', 'cum_int_n_half',
           'dis_to_hub_area', 'dis_to_hub_int', 'relative_r_area', 'relative_r_int']

for i in feature:
    f = i
    plt.subplots(figsize=(6, 4))
    weights1 = np.ones_like(data_sample[f]) / len(data_sample)
    weights2 = np.ones_like(data_WT[f]) / len(data_WT)
    weights3 = np.ones_like(data_control1[f]) / len(data_control1)
    plt.hist([data_sample[f], data_WT[f], data_control1[f]], weights=[weights1, weights2, weights3], color=hist_colors,
             edgecolor=(0.2, 0.2, 0.2), label=sample, bins=20)
    plt.xlabel(f)
    plt.ylabel('Probability')
    plt.legend()
    plt.savefig('%s/%s_%s.pdf' % (save_folder, f, sample))
    plt.close()

# single value feature-p
print("Plotting single value feature-p...")

for f in feature:
    dis.plot_minus_ln_p(5, 250, 50, f, data, control, sample, save_folder)

# multiple value feature
print("Plotting multiple value feature...")
data_m_sample = pd.DataFrame({'area_ind_ecDNA_m': data_sample_area_ind_ecDNA,
                              'area_ratio_ind_ecDNA_m': data_sample_area_ratio_ind_ecDNA,
                              'mean_int_ind_ecDNA_m': data_sample_mean_int_ind_ecDNA,
                              'sample': [sample]*len(data_sample_area_ind_ecDNA)})
data_m_control1 = pd.DataFrame({'area_ind_ecDNA_m': data_control1_area_ind_ecDNA,
                                'area_ratio_ind_ecDNA_m': data_control1_area_ratio_ind_ecDNA,
                                'mean_int_ind_ecDNA_m': data_control1_mean_int_ind_ecDNA,
                                'sample': [control[0]]*len(data_control1_area_ind_ecDNA)})
data_m_WT = pd.DataFrame({'area_ind_ecDNA_m': data_WT_area_ind_ecDNA,
                          'area_ratio_ind_ecDNA_m': data_WT_area_ratio_ind_ecDNA,
                          'mean_int_ind_ecDNA_m': data_WT_mean_int_ind_ecDNA,
                          'sample': ['WT']*len(data_WT_area_ind_ecDNA)})
data_m_control1_select = data_m_control1[data_m_control1['area_ind_ecDNA_m'] > 100].copy()
data_m_sample_select = data_m_sample[data_m_sample['area_ind_ecDNA_m'] > 100].copy()
data_m_WT_select = data_m_WT[data_m_WT['area_ind_ecDNA_m'] > 100].copy()
data_m = pd.concat([data_m_sample_select, data_m_WT_select, data_m_control1_select], axis=0, ignore_index=True)

feature = ['area_ind_ecDNA_m', 'mean_int_ind_ecDNA_m', 'area_ratio_ind_ecDNA_m']
for f in feature:
    plt.subplots(figsize=(6, 4))
    weights1 = np.ones_like(data_m_sample_select[f]) / len(data_m_sample_select)
    weights2 = np.ones_like(data_m_WT_select[f]) / len(data_m_WT_select)
    weights3 = np.ones_like(data_m_control1_select[f]) / len(data_m_control1_select)
    plt.hist([data_m_sample_select[f], data_m_WT_select[f], data_m_control1_select[f]], weights=[weights1, weights2, weights3],
             color=hist_colors, edgecolor=(0.2, 0.2, 0.2), label=sample, bins=20)
    plt.xlabel(f)
    plt.ylabel('Probability')
    plt.legend()
    plt.savefig('%s/%s_%s.pdf' % (save_folder, f, sample))
    plt.close()

# multiple value feature-p
print("Plotting multiple value feature-p...")
for f in feature:
    dis.plot_minus_ln_p(5, 250, 50, f, data_m, control, sample, save_folder)

# scatter
print("Plotting scatter images...")
feature_x_lst = ['mean_int_ecDNA', 'total_int_ecDNA', 'total_area_ecDNA', 'total_int_DNAFISH']
feature_y_lst = ['total_area_ratio_ecDNA', 'total_area_ecDNA']

sns.set_palette(sns.color_palette(line_colors))
data_part = data[(data['sample'] == sample) | (data['sample'] == 'WT')].copy()

for feature_x in feature_x_lst:
    for feature_y in feature_y_lst:
        ax1 = sns.jointplot(data=data_part, x=feature_x, y=feature_y, hue='sample', alpha=0.7, s=10)
        plt.savefig('%scomparison_of_%s_vs_%s.pdf' % (save_folder, feature_x, feature_y))
        plt.close()

# auto correlation curve
print("Plotting auto correlation curve...")
feature = ['g']
for f in feature:
    x = np.arange(0, 101, 1)
    x_label = 'r'

    number_nuclear1 = len(data_sample)
    number_nuclear2 = len(data_WT)

    mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data_sample[f].tolist())
    mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(data_WT[f].tolist())

    plt.subplots(figsize=(6, 4))
    for i in range(len(data_sample)):
        plt.plot(x, data_sample[f][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
    for i in range(len(data_WT)):
        plt.plot(x, data_WT[f][i], alpha=0.05, color=[line_colors[1][j]+0.05 for j in range(len(line_colors[1]))])
    plt.plot(x, mean_curve1, color=line_colors[0], label='%s, n=%s' % (sample, number_nuclear1))
    plt.plot(x, mean_curve2, color=line_colors[1], label='%s, n=%s' % ('WT', number_nuclear2))
    plt.plot(x, ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.xlabel(x_label)
    plt.ylabel(f)
    plt.legend()
    plt.savefig('%s/%s_%s.pdf' % (save_folder, f, sample))
    plt.close()

# radial curve
print("Plotting radial curve...")
feature = ['radial_curve_DNAFISH']
for f in feature:
    x = np.arange(0.01, 0.99, 0.01)
    x_label = 'relative r'

    number_nuclear1 = len(data_sample)
    number_nuclear2 = len(data_WT)

    mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data_sample[f].tolist())
    mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(data_WT[f].tolist())

    plt.subplots(figsize=(6, 4))
    for i in range(len(data_sample)):
        plt.plot(x, data_sample[f][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
    for i in range(len(data_WT)):
        plt.plot(x, data_WT[f][i], alpha=0.05, color=[line_colors[1][j]+0.05 for j in range(len(line_colors[1]))])
    plt.plot(x, mean_curve1, color=line_colors[0], label='%s, n=%s' % (sample, number_nuclear1))
    plt.plot(x, mean_curve2, color=line_colors[1], label='%s, n=%s' % ('WT', number_nuclear2))
    plt.plot(x, ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.xlabel(x_label)
    plt.ylabel(f)
    plt.legend()
    plt.savefig('%s/%s_%s.pdf' % (save_folder, f, sample))
    plt.close()

# angle curve
print("Plotting angle curve...")
feature = ['angle_curve_DNAFISH']
for f in feature:
    x = np.arange(0, 360, 1)
    x_label = 'degree'

    number_nuclear1 = len(data_sample)
    number_nuclear2 = len(data_WT)

    mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data_sample[f].tolist())
    mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(data_WT[f].tolist())

    plt.subplots(figsize=(6, 4))
    for i in range(len(data_sample)):
        plt.plot(x, data_sample[f][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
    for i in range(len(data_WT)):
        plt.plot(x, data_WT[f][i], alpha=0.05, color=[line_colors[1][j]+0.05 for j in range(len(line_colors[1]))])
    plt.plot(x, mean_curve1, color=line_colors[0], label='%s, n=%s' % (sample, number_nuclear1))
    plt.plot(x, mean_curve2, color=line_colors[1], label='%s, n=%s' % ('WT', number_nuclear2))
    plt.plot(x, ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.xlabel(x_label)
    plt.ylabel(f)
    plt.legend()
    plt.savefig('%s/%s_%s.pdf' % (save_folder, f, sample))
    plt.close()

# cumulative curve
print("Plotting cumulative curve...")
feature = ['cum_int_ind_ecDNA_filled', 'cum_area_ind_ecDNA_filled',
           'cum_area_ratio_ind_ecDNA_filled', 'percentage_area_curve_ecDNA',
           'percentage_int_curve_ecDNA', 'percentage_area_ratio_curve_ecDNA']
for f in feature:
    x_label = 'number of ecDNA hub'

    data_sample_filter = data_sample[data_sample['n_ecDNA'] > 5].copy().reset_index()
    data_control1_filter = data_control1[data_control1['n_ecDNA'] > 5].copy().reset_index()
    data_WT_filter = data_WT[data_WT['n_ecDNA'] > 5].copy().reset_index()

    number_nuclear1 = len(data_sample_filter)
    number_nuclear2 = len(data_WT_filter)
    number_nuclear3 = len(data_control1_filter)

    mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(dat.list_fill_with_last_num(data_sample_filter[f].tolist()))
    mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(dat.list_fill_with_last_num(data_WT_filter[f].tolist()))
    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(dat.list_fill_with_last_num(data_control1_filter[f].tolist()))

    plt.subplots(figsize=(6, 4))
    for i in range(len(data_sample)):
        plt.plot(data_sample[f][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
    for i in range(len(data_WT)):
        plt.plot(data_WT[f][i], alpha=0.05, color=[line_colors[1][j]+0.05 for j in range(len(line_colors[1]))])
    for i in range(len(data_control1)):
        plt.plot(data_control1[f][i], alpha=0.05, color=[line_colors[2][j]+0.05 for j in range(len(line_colors[2]))])
    plt.plot(mean_curve1, color=line_colors[0], label='%s, n=%s' % (sample, number_nuclear1))
    plt.plot(mean_curve2, color=line_colors[1], label='%s, n=%s' % ('WT', number_nuclear2))
    plt.plot(mean_curve3, color=line_colors[2], label='%s, n=%s' % ('Ctrl', number_nuclear3))
    plt.plot(ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(ci_lower3, color=line_colors[2], linestyle='--', linewidth=0.5)
    plt.plot(ci_higher3, color=line_colors[2], linestyle='--', linewidth=0.5)
    plt.xlabel(x_label)
    plt.ylabel(f)
    plt.legend()
    # plt.show()
    plt.savefig('%s/%s_%s.pdf' % (save_folder, f, sample))
    plt.close()

print("DONE!")

