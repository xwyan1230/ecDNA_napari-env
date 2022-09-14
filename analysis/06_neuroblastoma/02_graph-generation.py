import shared.segmentation as seg
from skimage.measure import regionprops, label
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, medial_axis
import shared.objects as obj
import shared.image as ima
from skimage import segmentation
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
group = 'ABC'
save_path = '%sv1_img/%s/' % (master_folder, group)
if not os.path.exists(save_path):
    os.makedirs(save_path)
version = 1
cmap = matplotlib.cm.get_cmap('Spectral')
compare_with = 'group'

data = pd.read_csv('%s%s_v%s_mean.txt' % (master_folder, group, version), na_values=['.'], sep='\t')

sample_lst = list(set(data[compare_with].tolist()))
sample_size = len(sample_lst)
rgba = [cmap(i) for i in np.arange(0, 1, 1/sample_size)]

feature_str = ['radial_curve_DNAFISH', 'radial_curve_nuclear', 'cum_area_ind_ecDNA',
               'cum_area_ratio_ind_ecDNA', 'percentage_area_curve_ecDNA', 'percentage_area_ratio_curve_ecDNA']
for f in feature_str:
    data[f] = [dat.str_to_float(data[f][i]) for i in range(len(data))]
data['radial_curve'] = [list(np.array(data['radial_curve_DNAFISH'].tolist()[i])/np.array(data['radial_curve_nuclear'].tolist()[i])) for i in range(len(data))]
data['area_ratio_per_ecDNA'] = data['total_area_ratio_ecDNA']/data['n_ecDNA']
# data = data[data['cum_area_n_half'] >= 6]

features = ['area_nuclear', 'n_ecDNA', 'total_area_ecDNA', 'total_area_ratio_ecDNA', 'max_area_ecDNA',
            'max_area_ratio_ecDNA', 'radial_center', 'radial_edge', 'relative_r_area', 'dis_to_hub_area',
            'percentage_area_ratio_n_half', 'percentage_area_n_half', 'cum_area_n_half', 'cum_area_ratio_n_half',
            'area_ratio_per_ecDNA']

data_A = data[data['group'] == 'A'].copy().reset_index(drop=True)
data_B = data[data['group'] == 'B'].copy().reset_index(drop=True)
data_C = data[data['group'] == 'C'].copy().reset_index(drop=True)

"""data = pd.DataFrame()
data = pd.concat([data_A, data_B, data_C], axis=0).reset_index(drop=True)"""
for f in features:
    histrange = [0, max(data[f].tolist())]
    weights1 = np.ones_like(data_A[f]) / len(data_A)
    weights2 = np.ones_like(data_B[f]) / len(data_B)
    weights3 = np.ones_like(data_C[f]) / len(data_C)
    plt.hist([data_A[f], data_B[f], data_C[f]], weights=[weights1, weights2, weights3], color=rgba,
             edgecolor=(0.2, 0.2, 0.2), label=['A', 'B', 'C'], bins=30)
    plt.xlabel(f)
    plt.ylabel('Probability')
    plt.legend()
    plt.savefig('%s%s_hist.pdf' % (save_path, f))
    plt.close()

    plt.subplots(figsize=(sample_size, 6))
    ax = sns.violinplot(x=compare_with, y=f, data=data)
    plt.savefig('%s%s.pdf' % (save_path, f))
    plt.close()

    sns.set_palette(sns.color_palette(rgba))
    sns.displot(data, x=f, hue="group", multiple="dodge")
    plt.savefig('%s%s_dodge.pdf' % (save_path, f))
    plt.close()

    sns.set_palette(sns.color_palette(rgba))
    sns.displot(data, x=f, hue="group", kind='kde', fill='True')
    plt.savefig('%s%s_kde.pdf' % (save_path, f))
    plt.close()

# radial curve
feature = ['radial_curve_DNAFISH', 'radial_curve_nuclear', 'radial_curve']
for f in feature:
    plt.subplots(figsize=(6, 4))
    x = np.arange(0.01, 0.99, 0.01)
    x_label = 'relative r'
    for sample in sample_lst:
        if len(data[data[compare_with] == sample]) != 0:
            data_sample = data[data[compare_with] == sample].copy().reset_index(drop=True)

            number_nuclear1 = len(data_sample)
            mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data_sample[f].tolist())
            plt.plot(x, mean_curve1, color=rgba[sample_lst.index(sample)], label='%s, n=%s' % (sample, number_nuclear1))
            # plt.plot(x, ci_lower1, color=rgba[sample_lst.index(sample)], linestyle='--', linewidth=0.5)
            # plt.plot(x, ci_higher1, color=rgba[sample_lst.index(sample)], linestyle='--', linewidth=0.5)

    plt.xlabel(x_label)
    plt.ylabel(f)
    plt.ylim([-1, 2])
    plt.legend()
    plt.savefig('%s/%s.pdf' % (save_path, f))
    plt.close()

# cumulative curve
feature = ['cum_area_ind_ecDNA', 'cum_area_ratio_ind_ecDNA', 'percentage_area_curve_ecDNA',
           'percentage_area_ratio_curve_ecDNA']
for f in feature:
    plt.subplots(figsize=(6, 4))
    x_label = 'number of ecDNA hub'

    for sample in sample_lst:
        if len(data[data[compare_with] == sample]) != 0:
            data_sample = data[data[compare_with] == sample].copy().reset_index()

            number_nuclear1 = len(data_sample)
            mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(dat.list_fill_with_last_num(data_sample[f].tolist()))
            plt.plot(mean_curve1, color=rgba[sample_lst.index(sample)], label='%s, n=%s' % (sample, number_nuclear1))
            # plt.plot(ci_lower1, color=rgba[sample_lst.index(sample)], linestyle='--', linewidth=0.5)
            # plt.plot(ci_higher1, color=rgba[sample_lst.index(sample)], linestyle='--', linewidth=0.5)

    plt.xlabel(x_label)
    plt.ylabel(f)
    if f in ['percentage_area_curve_ecDNA', 'percentage_area_ratio_curve_ecDNA']:
        plt.ylim([0, 1.2])
    plt.legend()
    plt.savefig('%s/%s.pdf' % (save_path, f))
    plt.close()

"""# radial curve
for sample in sample_lst:
    data_sample = data[data['sample'] == sample].copy().reset_index()
    feature = ['radial_curve_DNAFISH', 'radial_curve_nuclear', 'radial_curve']
    for f in feature:
        x = np.arange(0.01, 0.99, 0.01)
        x_label = 'relative r'

        number_nuclear1 = len(data_sample)
        mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data_sample[f].tolist())

        plt.subplots(figsize=(6, 4))
        for i in range(len(data_sample)):
            plt.plot(x, data_sample[f][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
        plt.plot(x, mean_curve1, color=line_colors[0], label='%s, n=%s' % (sample, number_nuclear1))
        plt.plot(x, ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
        plt.plot(x, ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
        plt.xlabel(x_label)
        plt.ylabel(f)
    plt.savefig('%s/%s_%s.pdf' % (save_path, sample, f))
    plt.close()

    # angle curve
    feature = ['angle_curve_DNAFISH']
    for f in feature:
        x = np.arange(0, 360, 1)
        x_label = 'degree'

        number_nuclear1 = len(data_sample)
        mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data_sample[f].tolist())

        plt.subplots(figsize=(6, 4))
        for i in range(len(data_sample)):
            plt.plot(x, data_sample[f][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
        plt.plot(x, mean_curve1, color=line_colors[0], label='%s, n=%s' % (sample, number_nuclear1))
        plt.plot(x, ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
        plt.plot(x, ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
        plt.xlabel(x_label)
        plt.ylabel(f)
        plt.savefig('%s/%s_%s.pdf' % (save_path, sample, f))
        plt.close()

    # cumulative curve
    feature = ['cum_area_ind_ecDNA_filled', 'cum_area_ratio_ind_ecDNA_filled', 'percentage_area_curve_ecDNA',
               'percentage_area_ratio_curve_ecDNA']
    for f in feature:
        x_label = 'number of ecDNA hub'

        number_nuclear1 = len(data_sample)
        mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(dat.list_fill_with_last_num(data_sample[f].tolist()))

        plt.subplots(figsize=(6, 4))
        for i in range(len(data_sample)):
            plt.plot(data_sample[f][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
        plt.plot(mean_curve1, color=line_colors[0], label='%s, n=%s' % (sample, number_nuclear1))
        plt.plot(ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
        plt.plot(ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
        plt.xlabel(x_label)
        plt.ylabel(f)
        plt.savefig('%s/%s_%s.pdf' % (save_path, sample, f))
        plt.close()"""

print("DONE!")
