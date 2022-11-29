import pandas as pd
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt
import matplotlib
import shared.math as mat
import seaborn as sns
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221017_periphery-localization_analysis/20220826_neuroblastoma/A/"
"""sample_lst = ['AU', 'BT', 'BZ', 'C', 'CF', 'CN', 'DK', 'DM', 'DY', 'DZ', 'EF', 'EH', 'ER', 'F', 'FI',
              'FL', 'FP', 'FV', 'FW', 'GA', 'GY', 'GZ', 'HC', 'HM', 'IP', 'IU', 'IX', 'JZ', 'KI', 'KJ',
              'KV', 'LM', 'ME', 'MU']  # C"""
# sample_lst = ['C', 'DZ', 'EF', 'F', 'FL', 'JZ', 'LM', 'ME']  # C select
sample_lst = ['AC', 'AK', 'CB', 'EB', 'FR', 'GL', 'GO', 'IC', 'NE', 'V']  # A
# sample_lst = ['EB', 'GO', 'NE']  # A select
group = 'A'
version = 1
fig_size = (12, 9)
# A (12, 9); C (24, 9)
filter = 'absolute'
save_folder = '%sfig/%s/' % (master_folder, filter)
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

# default setting
hist_colors = [(0.95, 0.50, 0.50), (0.90, 0.90, 0.90), (0.50, 0.90, 0.90)]
line_colors = [(0.85, 0.35, 0.25), (0.30, 0.30, 0.30), (0.30, 0.70, 0.70)]
line_color_4 = [(0.85, 0.35, 0.25), (0.30, 0.30, 0.30), (0.30, 0.70, 0.70), (0.95, 0.85, 0)]
cmap = matplotlib.cm.get_cmap('Spectral')
rgba = [cmap(i) for i in np.arange(0, 1, 1/len(sample_lst))]
sns.set_palette(sns.color_palette(rgba))
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
data['cluster_degree'] = 1/data['normalized_cluster']
data['radial_curve_normalized'] = \
    [list(np.array(data['radial_curve_DNAFISH'][i])/np.array(data['radial_curve_nuclear'][i])) for i in range(len(data))]

feature = 'radial_distance'
median_lst = []
for sample in sample_lst:
    temp = data[data['sample'] == sample][feature].median()
    median_lst.append(temp)
sample_lst_sorted = [x for _, x in sorted(zip(median_lst, sample_lst))]
print(sample_lst_sorted)
df = pd.DataFrame({'sample': sample_lst_sorted, 'radial_distance': sorted(median_lst)})
print(df)

data_sub = pd.DataFrame(columns=data.columns)
for sample in sample_lst_sorted:
    data_sample = data[data['sample'] == sample].copy()
    data_sub = pd.concat([data_sub, data_sample], ignore_index=True, axis=0)

# box plot
feature = ['radial_distance']
for f in feature:
    median_f = data_sub[f].median()
    plt.subplots(figsize=fig_size)
    ax = sns.boxplot(data=data_sub, y=f, x="sample", showmeans=True, showfliers=False)
    medians = []
    nobs = []
    for sample in sample_lst_sorted:
        medians.append(data_sub[data_sub['sample'] == sample][f].median())
        nobs.append(len(data_sub[data_sub['sample'] == sample]))
    nobs = [str(x) for x in nobs]
    nobs = ["n: " + i for i in nobs]
    pos = range(len(nobs))
    for tick, label in zip(pos, ax.get_xticklabels()):
        ax.text(pos[tick], medians[tick] + 0.1*median_f, nobs[tick], horizontalalignment='center', size=10, color='black',
                weight='semibold')
    plt.axhline(y=0, color='black', linestyle='--')
    plt.savefig('%s/%s-sample-%s.pdf' % (save_folder, group, f))
    plt.show()

feature = ['n_ecDNA', 'total_area_ratio_ecDNA']
for f in feature:
    median_f = data_sub[f].median()
    plt.subplots(figsize=fig_size)
    ax = sns.boxplot(data=data_sub, y=f, x="sample", showmeans=True, showfliers=False)
    medians = []
    nobs = []
    for sample in sample_lst_sorted:
        medians.append(data_sub[data_sub['sample'] == sample][f].median())
        nobs.append(len(data_sub[data_sub['sample'] == sample]))
    nobs = [str(x) for x in nobs]
    nobs = ["n: " + i for i in nobs]
    pos = range(len(nobs))
    for tick, label in zip(pos, ax.get_xticklabels()):
        ax.text(pos[tick], medians[tick] + 0.1*median_f, nobs[tick], horizontalalignment='center', size=10, color='black',
                weight='semibold')
    plt.savefig('%s/%s-sample-%s.pdf' % (save_folder, group, f))
    plt.show()

# radial curve
print("Plotting radial curve...")
x = np.arange(0.01, 0.99, 0.01)
x_label = 'relative r'
up_limit = 90
plt.subplots(figsize=(12, 9))
for sample in sample_lst_sorted:
    data_sample = data[data['sample'] == sample].copy()
    number_nuclear = len(data_sample)
    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data_sample['radial_curve_normalized'].tolist())
    plt.plot(x[:up_limit], mean_curve3[:up_limit], color=rgba[sample_lst_sorted.index(sample)], label='%s, n=%s' % (sample, number_nuclear))
    plt.plot(x[:up_limit], ci_lower3[:up_limit], color=rgba[sample_lst_sorted.index(sample)], linestyle='--', linewidth=0.5)
    plt.plot(x[:up_limit], ci_higher3[:up_limit], color=rgba[sample_lst_sorted.index(sample)], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.axhline(y=1, color='black', linestyle='--')
plt.ylim([0.75, 1.25])
plt.ylabel('radial_curve_normalized')
plt.legend()
plt.savefig('%s/%s-sample-radial_curve_normalized.pdf' % (save_folder, group))
plt.show()

# radial distribution vs cluster degree
data_filtered = data[data['n_ecDNA'] > 1].copy().reset_index(drop=True)
data_sub_filtered = data_sub[data_sub['n_ecDNA'] > 1].copy().reset_index(drop=True)

feature = ['cluster_degree', 'normalized_cluster']
for f in feature:
    median_f = data_sub[f].median()
    plt.subplots(figsize=fig_size)
    ax = sns.boxplot(data=data_sub_filtered, y=f, x="sample", showmeans=True, showfliers=False)
    medians = []
    nobs = []
    for sample in sample_lst_sorted:
        medians.append(data_sub_filtered[data_sub_filtered['sample'] == sample][f].median())
        nobs.append(len(data_sub_filtered[data_sub_filtered['sample'] == sample]))
    nobs = [str(x) for x in nobs]
    nobs = ["n: " + i for i in nobs]
    pos = range(len(nobs))
    for tick, label in zip(pos, ax.get_xticklabels()):
        ax.text(pos[tick], medians[tick] + 0.1*median_f, nobs[tick], horizontalalignment='center', size=10, color='black',
                weight='semibold')
    plt.savefig('%s/%s-sample-%s.pdf' % (save_folder, group, f))
    plt.show()

plt.subplots(figsize=(12, 9))
for sample in sample_lst_sorted:
    data_sample = data_filtered[data_filtered['sample'] == sample].copy()
    number_nuclear = len(data_sample)
    plt.scatter(data_sample['total_area_ecDNA']**0.5, data_sample['dis_to_hub_area_v2'],
                c=rgba[sample_lst_sorted.index(sample)], label='%s, n=%s' % (sample, number_nuclear))
plt.xlabel('total_area_ecDNA_sqrt')
plt.ylabel('dis_to_hub_area_v2')
plt.legend()
plt.savefig('%s/%s-sample-dis_to_hub_area_v2-vs-total_area_ecDNA_sqrt.pdf' % (save_folder, group))
plt.show()

k_cluster = []
for sample in sample_lst_sorted:
    plt.subplots(figsize=(12, 9))
    data_sample = data_filtered[data_filtered['sample'] == sample].copy()
    number_nuclear = len(data_sample)
    plt.scatter(data_sample['total_area_ecDNA'] ** 0.5, data_sample['dis_to_hub_area_v2'],
                c=rgba[sample_lst_sorted.index(sample)], label='%s, n=%s' % (sample, number_nuclear))
    _, int_fit_r2_sample, int_fit_a_sample = \
        mat.fitting_linear_b0(list(data_sample['total_area_ecDNA'] ** 0.5), data_sample['dis_to_hub_area_v2'].tolist())
    k_cluster.append(int_fit_a_sample[0])
    x = np.arange(0, 1, 0.01)
    y_sample = int_fit_a_sample * x
    plt.plot(x, y_sample, color='black', linestyle='--', label='r2=%s, a=%s' % (int_fit_r2_sample, int_fit_a_sample[0]))
    plt.xlabel('total_area_ecDNA_sqrt')
    plt.ylabel('dis_to_hub_area_v2')
    plt.legend()
    plt.savefig('%s/%s-sample_%s-dis_to_hub_area_v2-vs-total_area_ratio_ecDNA_sqrt.pdf' % (save_folder, group, sample))
    plt.show()

df['k_cluster'] = k_cluster
plt.subplots(figsize=(12, 9))
plt.scatter(df['radial_distance'], df['k_cluster'], color='red', label='group %s, n=%s' % (group, len(sample_lst)))
_, int_fit_r2_sample, int_fit_a_sample, int_fit_b_sample = \
    mat.fitting_linear(df['radial_distance'].tolist(), df['k_cluster'].tolist())
x = np.arange(min(df['radial_distance'])-0.01, max(df['radial_distance'])+0.01, 0.01)
y_sample = int_fit_a_sample * x + int_fit_b_sample
plt.plot(x, y_sample, color='black', linestyle='--', label='r2=%s' % int_fit_r2_sample)
plt.xlabel('radial_distance')
plt.ylabel('k_cluster')
plt.legend()
plt.savefig('%s/%s-sample-radial_distance-vs-k_cluster.pdf' % (save_folder, group))
plt.show()
df.to_csv('%s%s_df_v%s.txt' % (save_folder, group, version), index=False, sep='\t')

plt.subplots(figsize=(12, 9))
for sample in sample_lst_sorted:
    data_sample = data_filtered[data_filtered['sample'] == sample].copy()
    number_nuclear = len(data_sample)
    plt.scatter(data_sample['radial_distance'], data_sample['normalized_cluster'],
                color=rgba[sample_lst_sorted.index(sample)], label='%s, n=%s' % (sample, number_nuclear))
plt.xlabel('radial_distance')
plt.ylabel('cluster_degree')
plt.legend()
plt.savefig('%s/%s-sample-radial_distance-vs-normalized_cluster.pdf' % (save_folder, group))
plt.show()

print("DONE!")