import skimage.io as skio
import pandas as pd
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.math as mat
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import shared.dataframe as dat
import os
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230119_model_clustering/"
data_dir = "%stxt/dataset5/60/" % master_folder
output_dir = "%sfigures1/" % master_folder

samples = ['0_5', '1_5', '2_5', '5_5', '10_5', '25_5', '50_5', '75_5', '100_5', '200_5', '300_5', '400_5', '500_5', '1000_5',
           '2000_5', '3000_5', '4000_5', '5000_5']

cmap = mpl.cm.get_cmap('Spectral')
x = np.arange(0, 1, 1/len(samples))
line_color = [cmap(i) for i in x]

data = pd.DataFrame()
for i in samples:
    df = pd.read_csv(("%s%s_cluster.txt" % (data_dir, i)), na_values=['.'], sep='\t')
    feature = ['cum_area_ind_ecDNA_filled', 'percentage_area_curve_ecDNA', 'g']
    for f in feature:
        df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]
    data = pd.concat([data, df], axis=0)

data['r'] = np.sqrt(data['area_nuclear']/math.pi)
data['total_area_ecDNA_sqrt'] = np.sqrt(data['total_area_ecDNA']/math.pi)
data['total_area_ecDNA_sqrt_normalized'] = data['total_area_ecDNA_sqrt']/data['r']
data['dis_to_hub_area_v2_normalized'] = data['dis_to_hub_area_v2']/data['r']
data['cum_area_ind_ecDNA_filled_all'] = dat.list_fill_with_last_num(data['cum_area_ind_ecDNA_filled'].tolist())
data_filter = data[(data['copy_num'] < 50) & (data['total_area_ecDNA_sqrt_normalized'] > 0.1)].copy().reset_index(drop=True)
# data_filter = data[data['copy_num'] < 50].copy().reset_index(drop=True)
data = data_filter

low_percentage = []  # <0.2
mid_percentage = []  # [0.2, 0.8)
high_percentage = []  # >=0.8

for i in range(len(samples)):
    coefficient = int(samples[i].split('_')[0])
    data_sample = data[data['coefficient'] == coefficient].copy().reset_index(drop=True)
    low_percentage.append(len(data_sample[data_sample['dis_to_hub_area_v2_normalized'] < 0.2])*1.0/len(data_sample))
    mid_percentage.append(len(data_sample[(data_sample['dis_to_hub_area_v2_normalized'] >= 0.2) & (data_sample['dis_to_hub_area_v2_normalized'] < 0.8)]) * 1.0 / len(data_sample))
    high_percentage.append(len(data_sample[data_sample['dis_to_hub_area_v2_normalized'] >= 0.8])*1.0/len(data_sample))
plt.subplots(figsize=(12, 9))
plt.bar(samples, low_percentage, color=line_color)
plt.xlabel('sample coefficient')
plt.ylabel('percentage of low dis_to_hub_area_v2_normalized')
plt.savefig('%slow_percentage.pdf' % output_dir)
plt.close()

plt.subplots(figsize=(12, 9))
plt.bar(samples, mid_percentage, color=line_color)
plt.xlabel('sample coefficient')
plt.ylabel('percentage of mid dis_to_hub_area_v2_normalized')
plt.savefig('%smid_percentage.pdf' % output_dir)
plt.close()

plt.subplots(figsize=(12, 9))
plt.bar(samples, high_percentage, color=line_color)
plt.xlabel('sample coefficient')
plt.ylabel('percentage of high dis_to_hub_area_v2_normalized')
plt.savefig('%shigh_percentage.pdf' % output_dir)
plt.close()

d_range = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
mean_df = []
plt.subplots(figsize=(12, 9))
for i in range(len(samples)):
    coefficient = int(samples[i].split('_')[0])
    data_sample = data[data['coefficient'] == coefficient].copy().reset_index(drop=True)
    temp = []
    for j in range(len(d_range)-1):
        temp.append(len(data_sample[(data_sample['dis_to_hub_area_v2_normalized'] >= d_range[j]) & (data_sample['dis_to_hub_area_v2_normalized'] < d_range[j+1])]) * 1.0 / len(data_sample))
    mean_df.append(temp)
    plt.plot(d_range[1:], temp, color=line_color[i], label='%s, n=%s' % (coefficient, len(data_sample)))
plt.xlabel('dis_to_hub_area_v2_normalized')
plt.ylabel('dis_curve')
plt.legend(loc=2)
plt.savefig('%sdis_curve.pdf' % output_dir)
plt.show()

# heatmap
df = pd.DataFrame(mean_df)
df.columns = ['%.1f' % elem for elem in d_range][1:]
df.index = samples
plt.subplots(figsize=(12, 9))
sns.heatmap(df, cbar=0, linewidths=2, square=True, cmap='coolwarm')
plt.savefig('%sdis_heatmap.pdf' % output_dir)
plt.show()

dis = pd.DataFrame()
dis['coefficient'] = [int(x.split('_')[0]) for x in samples] * 6
dis['dis_to_hub'] = [0.1] * len(samples) + [0.3] * len(samples) + [0.5] * len(samples) + [0.7] * len(samples) + [0.9] * len(samples) + [1.1] * len(samples)
dis['percentage'] = df['0.2'].tolist() + df['0.4'].tolist() + df['0.6'].tolist() + df['0.8'].tolist() + df['1.0'].tolist() + df['1.2'].tolist()
print(dis.head())

plt.subplots(figsize=(20, 9))
sns.barplot(data=dis, x='coefficient', y='percentage', hue='dis_to_hub')
plt.savefig('%sdis_barplot.pdf' % output_dir)
plt.show()
