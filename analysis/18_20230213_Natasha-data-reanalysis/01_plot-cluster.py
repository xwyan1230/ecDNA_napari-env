import skimage.io as skio
import pandas as pd
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.math as mat
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import shared.dataframe as dat
import math
import os
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230213_analysis_Natasha-data-reanalysis/"
data_dir = "%stxt/20220708_Natasha_ColoDM_interphase/" % master_folder
output_dir = "%sfigures/20220708_Natasha_ColoDM_interphase/" % master_folder

samples = ['DMSO1', '15hrJQ1', 'DMSO2']
line_color = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25),  (0.30, 0.70, 0.70)]
sns.set_palette(sns.color_palette(line_color))

data = pd.DataFrame()
for i in samples:
    df = pd.read_csv(("%s%s.txt" % (data_dir, i)), na_values=['.'], sep='\t')
    feature = ['cum_area_ind_ecDNA_filled', 'percentage_area_curve_ecDNA', 'g']
    df['sample'] = [i] * len(df)
    for f in feature:
        df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]
    data = pd.concat([data, df], axis=0)

data['r'] = np.sqrt(data['area_nuclear']/math.pi)
data['total_area_ecDNA_sqrt'] = np.sqrt(data['total_area_ecDNA']/math.pi)
data['total_area_ecDNA_sqrt_normalized'] = data['total_area_ecDNA_sqrt']/data['r']
data['dis_to_hub_area_v2_normalized'] = data['dis_to_hub_area']/data['r']
data['cum_area_ind_ecDNA_filled_all'] = dat.list_fill_with_last_num(data['cum_area_ind_ecDNA_filled'].tolist())

# data_filter = data[data['dis_to_hub_area_v2_normalized'] < 1.2].copy().reset_index(drop=True)
data_filter = data[(data['dis_to_hub_area_v2_normalized'] < 1.2) & (data['total_area_ecDNA_sqrt_normalized'] > 0.3)].copy().reset_index(drop=True)
data = data_filter

plt.subplots(figsize=(4, 9))
sns.violinplot(data=data, x='sample', y='dis_to_hub_area_v2_normalized')
plt.savefig('%sdis_to_hub_area_v2_normalized.pdf' % output_dir)
plt.show()

# plot cluster for all the samples
plt.subplots(figsize=(12, 9))
for i in range(len(samples)):
    data_sample = data[data['sample'] == samples[i]].copy().reset_index(drop=True)
    sns.scatterplot(data=data_sample, x='total_area_ecDNA_sqrt_normalized', y='dis_to_hub_area_v2_normalized', alpha=1, color=line_color[i],
                    label=samples[i])
plt.xlim([0, 0.7])
plt.ylim([-0.05, 1.2]) # 1.05 for cen_r 0, 1.2 for cen_r>0
plt.legend()
plt.savefig('%sdis_to_hub_area_normralized_vs_total_area_sqrt_normalized.pdf' % output_dir)
plt.show()

# plot cluster for individual sample
for i in range(len(samples)):
    data_sample = data[data['sample'] == samples[i]].copy().reset_index(drop=True)
    plt.subplots(figsize=(12, 9))
    sns.scatterplot(data=data_sample, x='total_area_ecDNA_sqrt_normalized', y='dis_to_hub_area_v2_normalized', color=line_color[i])
    plt.xlim([0, 0.7])
    plt.ylim([-0.05, 1.2])
    plt.savefig('%s%s_dis_to_hub_area_normalized_vs_total_area_sqrt_normalized.pdf' % (output_dir, samples[i]))
    plt.close()

# plot g curve for all samples
x = np.arange(0, 76, 1)
x_label = 'r'
limit = 75
plt.subplots(figsize=(12, 9))

for i in range(len(samples)):
    data_sample = data[data['sample'] == samples[i]].copy().reset_index(drop=True)
    number_nuclear = len(data_sample)
    mean_curve, ci_lower, ci_higher = dat.mean_list(data_sample['g'].tolist())
    plt.plot(x[1:limit], mean_curve[1:limit], color=line_color[i], label='%s, n=%s' % (samples[i], number_nuclear))

plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylabel('g')
plt.legend()
plt.savefig('%sg.pdf' % output_dir)
plt.show()

# plot g curve for individual samples
for i in range(len(samples)):
    data_sample = data[data['sample'] == samples[i]].copy().reset_index(drop=True)
    number_nuclear = len(data_sample)
    mean_curve, ci_lower, ci_higher = dat.mean_list(data_sample['g'].tolist())
    plt.subplots(figsize=(12, 9))
    for j in range(len(data_sample)):
        plt.plot(x[1:limit], data_sample['g'][j][1:limit], alpha=0.02, color=line_color[i])
    plt.plot(x[1:limit], mean_curve[1:limit], color=line_color[i], label='%s, n=%s' % (samples[i], number_nuclear))
    plt.plot(x[1:limit], ci_lower[1:limit], color=line_color[i], linestyle='--', linewidth=0.5)
    plt.plot(x[1:limit], ci_higher[1:limit], color=line_color[i], linestyle='--', linewidth=0.5)

    plt.axhline(y=1, color='black', linestyle='--')
    plt.xlabel(x_label)
    plt.ylabel('g')
    plt.legend()
    plt.savefig('%s%s_g.pdf' % (output_dir, samples[i]))
    plt.close()

# plot percentage/cumulative curve for all samples
feature = ['cum_area_ind_ecDNA_filled_all', 'percentage_area_curve_ecDNA']
for f in feature:
    plt.subplots(figsize=(12, 9))
    x_label = 'number of ecDNA hub'
    for i in range(len(samples)):
        data_sample = data[data['sample'] == samples[i]].copy().reset_index(drop=True)
        for j in range(len(data_sample)):
            plt.plot(data_sample[f][j], alpha=0.01, color=line_color[i])
    for i in range(len(samples)):
        data_sample = data[data['sample'] == samples[i]].copy().reset_index(drop=True)
        number_nuclear = len(data_sample)
        mean_curve, ci_lower, ci_higher = dat.mean_list(dat.list_fill_with_last_num(data_sample[f].tolist()))
        plt.plot(mean_curve, color=line_color[i], label='%s, n=%s' % (samples[i], number_nuclear))
    plt.xlabel(x_label)
    plt.ylabel(f)
    plt.legend()
    plt.savefig('%s%s.pdf' % (output_dir, f))
    plt.show()

# plot percentage/cumulative curve for individual sample
feature = ['cum_area_ind_ecDNA_filled', 'percentage_area_curve_ecDNA']
for i in range(len(samples)):
    data_sample = data[data['sample'] == samples[i]].copy().reset_index(drop=True)
    number_nuclear = len(data_sample)
    for f in feature:
        plt.subplots(figsize=(12, 9))
        for j in range(len(data_sample)):
            plt.plot(data_sample[f][j], alpha=0.05, color=line_color[i])
        mean_curve, ci_lower, ci_higher = dat.mean_list(dat.list_fill_with_last_num(data_sample[f].tolist()))
        plt.plot(mean_curve, color=line_color[i], label='%s, n=%s' % (samples[i], number_nuclear))
        plt.plot(ci_lower, color=line_color[i], linestyle='--', linewidth=0.5)
        plt.plot(ci_higher, color=line_color[i], linestyle='--', linewidth=0.5)
        plt.xlabel(x_label)
        plt.ylabel(f)
        plt.legend()
        plt.savefig('%s%s_%s.pdf' % (output_dir, samples[i], f))
        plt.close()


d_range = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
mean_df = []
plt.subplots(figsize=(12, 9))
for i in range(len(samples)):
    data_sample = data[data['sample'] == samples[i]].copy().reset_index(drop=True)
    temp = []
    for j in range(len(d_range)-1):
        temp.append(len(data_sample[(data_sample['dis_to_hub_area_v2_normalized'] >= d_range[j]) & (data_sample['dis_to_hub_area_v2_normalized'] < d_range[j+1])]) * 1.0 / len(data_sample))
    mean_df.append(temp)

df_temp = pd.DataFrame(mean_df)
df_temp.columns = ['%.1f' % elem for elem in d_range][1:]
df_temp.index = samples

dis = pd.DataFrame()
dis['sample'] = list(samples) * 6
dis['dis_to_hub'] = [0.1] * len(samples) + [0.3] * len(samples) + [0.5] * len(samples) + [0.7] * len(samples) + [0.9] * len(samples) + [1.1] * len(samples)
dis['percentage'] = df_temp['0.2'].tolist() + df_temp['0.4'].tolist() + df_temp['0.6'].tolist() + df_temp['0.8'].tolist() + df_temp['1.0'].tolist() + df_temp['1.2'].tolist()
print(dis.head())

plt.subplots(figsize=(9, 9))
sns.barplot(data=dis, x='sample', y='percentage', hue='dis_to_hub')
plt.savefig('%sdis_barplot.pdf' % output_dir)
plt.show()