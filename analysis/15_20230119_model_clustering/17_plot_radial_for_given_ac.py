import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import matplotlib as mpl
import numpy as np
import math
from matplotlib import cm
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230119_model_clustering/"
data_dir = "%stxt/dataset1/" % master_folder
output_dir = "%sfigures2_d1/" % master_folder

sample = '200_5'
# interval = 5
# seps = np.arange(50, 101, interval)
# sep_feature = 'r'
interval = 25
seps = np.arange(0, 201, interval)
sep_feature = 'copy_num'

cmap = mpl.cm.get_cmap('Spectral')
x = np.arange(0, 1, 1/len(seps))
line_color = [cmap(i) for i in x]

data = pd.read_csv(("%s%s_cluster.txt" % (data_dir, sample)), na_values=['.'], sep='\t')
feature = ['radial_curve_DNAFISH']
for f in feature:
    data[f] = [dat.str_to_float(data[f][i]) for i in range(len(data))]
data['r'] = np.sqrt(data['area_nuclear']/math.pi)

# radial curve
x = np.arange(0.025, 1, 0.05)

mean_df = []
plt.subplots(figsize=(12, 9))
for i in range(len(seps)-1):
    coefficient = int(sample.split('_')[0])
    data_sample = data[(data[sep_feature] > seps[i]) & (data[sep_feature] <= seps[i+1])].copy().reset_index(drop=True)
    for j in range(len(data_sample)):
        plt.plot(x, data_sample['radial_curve_DNAFISH'][j], alpha=0.01, color=line_color[i])
    mean_curve, ci_lower, ci_higher = dat.mean_list(data_sample['radial_curve_DNAFISH'].tolist())
    mean_df.append(mean_curve)
    plt.plot(x, mean_curve, color=line_color[i], label='%s, n=%s' % (seps[i+1], len(data_sample)))
plt.xlabel('relative_r')
plt.ylabel('radial_curve')
plt.ylim([0, 30])
plt.legend(loc=2)
plt.savefig('%s%s_radial_curve.pdf' % (output_dir, sample))
plt.show()

# heatmap
df = pd.DataFrame(mean_df)
df.columns = ['%.3f' % elem for elem in x]
df.index = seps[1:]
plt.subplots(figsize=(12, 9))
sns.heatmap(df, cbar=0, linewidths=2, square=True, cmap='coolwarm')
plt.savefig('%s%s_radial_heatmap.pdf' % (output_dir, sample))
plt.show()