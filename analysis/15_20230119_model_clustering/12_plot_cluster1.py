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
data_dir = "%stxt/" % master_folder
output_dir = "%sfigures2/" % master_folder

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

"""# plot cluster for all the samples
plt.subplots(figsize=(12, 9))
for i in range(len(samples)):
    coefficient = int(samples[i].split('_')[0])
    data_sample = data[data['coefficient'] == coefficient].copy().reset_index(drop=True)
    sns.scatterplot(data=data_sample, x='total_area_ecDNA_sqrt_normalized', y='dis_to_hub_area_v2_normalized', alpha=0.5, color=line_color[i],
                    label=coefficient)
plt.xlim([-0.01, 0.5])
plt.ylim([-0.01, 1])
plt.legend()
plt.savefig('%sdis_to_hub_area_vs_total_area_sqrt.pdf' % output_dir)
plt.show()

# plot cluster for individual sample
for i in range(len(samples)):
    coefficient = int(samples[i].split('_')[0])
    data_sample = data[data['coefficient'] == coefficient].copy().reset_index(drop=True)
    plt.subplots(figsize=(12, 9))
    sns.scatterplot(data=data_sample, x='total_area_ecDNA_sqrt_normalized', y='dis_to_hub_area_v2_normalized', c=data_sample['copy_num'])
    plt.xlim([-0.01, 0.5])
    plt.ylim([-0.01, 1])
    plt.savefig('%s%s_dis_to_hub_area_vs_total_area_sqrt.pdf' % (output_dir, coefficient))
    plt.close()"""

plt.subplots(figsize=(12, 9))
sns.violinplot(data=data, x='coefficient', y='dis_to_hub_area_v2_normalized')
plt.savefig('%sdis_to_hub_area_v2_normalized.pdf' % output_dir)
plt.show()