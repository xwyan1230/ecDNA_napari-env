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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230119_model_clustering/"
data_dir = "%stxt/dataset6_dataset3-radial/" % master_folder
output_dir = "%sfigures3/" % master_folder

# samples = ['0-50_5', '10-50_5', '20-50_5', '30-50_5', '40-50_5', '50-50_5', '60-50_5', '70-50_5']
samples = ['40-0_5', '40-1_5', '40-2_5', '40-5_5', '40-10_5', '40-25_5', '40-50_5', '40-75_5', '40-100_5', '40-200_5', '40-300_5', '40-400_5', '40-500_5', '40-1000_5',
           '40-2000_5', '40-3000_5', '40-4000_5', '40-5000_5']
hue_order = samples

cmap = mpl.cm.get_cmap('Spectral')
x = np.arange(0, 1, 1/len(samples))
line_color = [cmap(i) for i in x]

data = pd.DataFrame()
for i in samples:
    cen_r = i.split('-')[0]
    sample = i.split('-')[1]
    coefficient = sample.split('_')[0]
    df = pd.read_csv(("%s%s/%s_radial.txt" % (data_dir, cen_r, sample)), na_values=['.'], sep='\t')
    feature = ['cum_area_ind_ecDNA_filled', 'percentage_area_curve_ecDNA', 'g']
    for f in feature:
        df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]
    data = pd.concat([data, df], axis=0)

data['r'] = np.sqrt(data['area_nuclear']/math.pi)
data['total_area_ecDNA_sqrt'] = np.sqrt(data['total_area_ecDNA']/math.pi)
data['total_area_ecDNA_sqrt_normalized'] = data['total_area_ecDNA_sqrt']/data['r']
data['dis_to_hub_area_v2_normalized'] = data['dis_to_hub_area_v2']/data['r']
data['cum_area_ind_ecDNA_filled_all'] = dat.list_fill_with_last_num(data['cum_area_ind_ecDNA_filled'].tolist())
df = data

# plot cluster for all the samples
plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    cen_r = int(hue_order[k].split('-')[0])
    sample = hue_order[k].split('-')[1]
    coefficient = int(sample.split('_')[0])
    data = df[(df['cen_r'] == cen_r) & (df['coefficient'] == coefficient)].copy().reset_index(drop=True)
    sns.scatterplot(data=data, x='total_area_ecDNA_sqrt_normalized', y='dis_to_hub_area_v2_normalized', alpha=0.5, color=line_color[k],
                    label=hue_order[k])
plt.xlim([0, 0.45])
plt.ylim([-0.05, 1.3]) # 1.05 for cen_r 0, 1.2 for cen_r>0
plt.legend()
plt.savefig('%sdis_to_hub_area_normralized_vs_total_area_sqrt_normalized.pdf' % output_dir)
plt.show()

# plot cluster for individual sample
for k in range(len(hue_order)):
    cen_r = int(hue_order[k].split('-')[0])
    sample = hue_order[k].split('-')[1]
    coefficient = int(sample.split('_')[0])
    data = df[(df['cen_r'] == cen_r) & (df['coefficient'] == coefficient)].copy().reset_index(drop=True)
    plt.subplots(figsize=(12, 9))
    sns.scatterplot(data=data, x='total_area_ecDNA_sqrt_normalized', y='dis_to_hub_area_v2_normalized', color=line_color[k])
    plt.xlim([0, 0.45])
    plt.ylim([-0.05, 1.3])
    plt.savefig('%s%s_dis_to_hub_area_normalized_vs_total_area_sqrt_normalized.pdf' % (output_dir, hue_order[k]))
    plt.close()

# plot g curve for all samples
x = np.arange(0, 76, 1)
x_label = 'r'
limit = 75
plt.subplots(figsize=(12, 9))

for k in range(len(hue_order)):
    cen_r = int(hue_order[k].split('-')[0])
    sample = hue_order[k].split('-')[1]
    coefficient = int(sample.split('_')[0])
    data = df[(df['cen_r'] == cen_r) & (df['coefficient'] == coefficient)].copy().reset_index(drop=True)
    number_nuclear = len(data)
    mean_curve, ci_lower, ci_higher = dat.mean_list(data['g'].tolist())
    plt.plot(x[1:limit], mean_curve[1:limit], color=line_color[k], label='%s, n=%s' % (hue_order[k], number_nuclear))

plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylabel('g')
plt.legend()
plt.savefig('%sg.pdf' % output_dir)
plt.show()

# plot percentage/cumulative curve for all samples
feature = ['cum_area_ind_ecDNA_filled_all', 'percentage_area_curve_ecDNA']
for f in feature:
    plt.subplots(figsize=(12, 9))
    x_label = 'number of ecDNA hub'
    for k in range(len(hue_order)):
        cen_r = int(hue_order[k].split('-')[0])
        sample = hue_order[k].split('-')[1]
        coefficient = int(sample.split('_')[0])
        data = df[(df['cen_r'] == cen_r) & (df['coefficient'] == coefficient)].copy().reset_index(drop=True)
        for j in range(len(data)):
            plt.plot(data[f][j], alpha=0.01, color=line_color[k])
    for k in range(len(hue_order)):
        cen_r = int(hue_order[k].split('-')[0])
        sample = hue_order[k].split('-')[1]
        coefficient = int(sample.split('_')[0])
        data = df[(df['cen_r'] == cen_r) & (df['coefficient'] == coefficient)].copy().reset_index(drop=True)
        number_nuclear = len(data)
        mean_curve, ci_lower, ci_higher = dat.mean_list(dat.list_fill_with_last_num(data[f].tolist()))
        plt.plot(mean_curve, color=line_color[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.xlabel(x_label)
    plt.ylabel(f)
    plt.legend()
    plt.savefig('%s%s.pdf' % (output_dir, f))
    plt.show()