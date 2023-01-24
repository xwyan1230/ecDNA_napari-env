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
import os
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230119_model_clustering/"
data_dir = "%stxt/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = '5000_5'
interval = 25
samples = np.arange(0, 201, interval)
print(samples)

cmap = mpl.cm.get_cmap('Spectral')
x = np.arange(0, 1, 1/len(samples))
line_color = [cmap(i) for i in x]

data = pd.read_csv(("%s%s_cluster.txt" % (data_dir, sample)), na_values=['.'], sep='\t')
feature = ['cum_area_ind_ecDNA_filled', 'percentage_area_curve_ecDNA', 'g']
for f in feature:
    data[f] = [dat.str_to_float(data[f][i]) for i in range(len(data))]
data['total_area_ecDNA_sqrt'] = np.sqrt(data['total_area_ecDNA'])
data['cum_area_ind_ecDNA_filled_all'] = dat.list_fill_with_last_num(data['cum_area_ind_ecDNA_filled'].tolist())

# plot cluster for all the samples
plt.subplots(figsize=(12, 9))
for i in range(len(samples)-1):
    data_sample = data[(data['copy_num'] > samples[i]) & (data['copy_num'] <= samples[i+1])].copy().reset_index(drop=True)
    sns.scatterplot(data=data_sample, x='total_area_ecDNA_sqrt', y='dis_to_hub_area_v2', alpha=0.5, color=line_color[i],
                    label=samples[i+1])
plt.xlim([10, 50])
plt.ylim([-5, 80])
plt.legend()
plt.savefig('%scopy_num_%s_dis_to_hub_area_vs_total_area_sqrt.pdf' % (output_dir, sample))
plt.show()

# plot g curve for all samples
x = np.arange(0, 76, 1)
x_label = 'r'
limit = 75
plt.subplots(figsize=(12, 9))

for i in range(len(samples)-1):
    data_sample = data[(data['copy_num'] > samples[i]) & (data['copy_num'] <= samples[i+1])].copy().reset_index(drop=True)
    number_nuclear = len(data_sample)
    mean_curve, ci_lower, ci_higher = dat.mean_list(data_sample['g'].tolist())
    plt.plot(x[1:limit], mean_curve[1:limit], color=line_color[i], label='%s, n=%s' % (samples[i+1], number_nuclear))

plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylabel('g')
plt.legend()
plt.savefig('%scopy_num_%s_g.pdf' % (output_dir, sample))
plt.show()

# plot percentage/cumulative curve for all samples
feature = ['cum_area_ind_ecDNA_filled_all', 'percentage_area_curve_ecDNA']
for f in feature:
    plt.subplots(figsize=(12, 9))
    x_label = 'number of ecDNA hub'
    for i in range(len(samples) - 1):
        data_sample = data[(data['copy_num'] > samples[i]) & (data['copy_num'] <= samples[i + 1])].copy().reset_index(
            drop=True)
        for j in range(len(data_sample)):
            plt.plot(data_sample[f][j], alpha=0.01, color=line_color[i])
    for i in range(len(samples) - 1):
        data_sample = data[(data['copy_num'] > samples[i]) & (data['copy_num'] <= samples[i + 1])].copy().reset_index(
            drop=True)
        number_nuclear = len(data_sample)
        mean_curve, ci_lower, ci_higher = dat.mean_list(dat.list_fill_with_last_num(data_sample[f].tolist()))
        plt.plot(mean_curve, color=line_color[i], label='%s, n=%s' % (samples[i+1], number_nuclear))
    plt.xlabel(x_label)
    plt.ylabel(f)
    plt.legend()
    plt.savefig('%scopy_num_%s_%s.pdf' % (output_dir, sample, f))
    plt.show()