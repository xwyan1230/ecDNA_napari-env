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
data_dir = "%stxt/dataset2/" % master_folder
output_dir = "%sfigures2_d2/" % master_folder

samples = ['0_5', '1_5', '2_5', '5_5', '10_5', '25_5', '50_5', '75_5', '100_5', '200_5', '300_5', '400_5', '500_5', '1000_5',
           '2000_5', '3000_5', '4000_5', '5000_5']

cmap = mpl.cm.get_cmap('Spectral')
x = np.arange(0, 1, 1/len(samples))
line_color = [cmap(i) for i in x]

data = pd.DataFrame()
for i in samples:
    df = pd.read_csv(("%s%s_cluster.txt" % (data_dir, i)), na_values=['.'], sep='\t')
    feature = ['radial_curve_DNAFISH']
    for f in feature:
        df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]
    data = pd.concat([data, df], axis=0)

data['r'] = np.sqrt(data['area_nuclear']/math.pi)

# radial curve
x = np.arange(0.025, 1, 0.05)

mean_df = []
plt.subplots(figsize=(12, 9))
for i in range(len(samples)):
    coefficient = int(samples[i].split('_')[0])
    data_sample = data[data['coefficient'] == coefficient].copy().reset_index(drop=True)
    for j in range(len(data_sample)):
        plt.plot(x, data_sample['radial_curve_DNAFISH'][j], alpha=0.01, color=line_color[i])
    mean_curve, ci_lower, ci_higher = dat.mean_list(data_sample['radial_curve_DNAFISH'].tolist())
    mean_df.append(mean_curve)
    plt.plot(x, mean_curve, color=line_color[i], label='%s, n=%s' % (coefficient, len(data_sample)))
plt.xlabel('relative_r')
plt.ylabel('radial_curve')
plt.ylim([0, 20])
plt.legend(loc=2)
plt.savefig('%sradial_curve.pdf' % output_dir)
plt.show()

# heatmap
df = pd.DataFrame(mean_df)
df.columns = ['%.3f' % elem for elem in x]
df.index = samples
plt.subplots(figsize=(12, 9))
sns.heatmap(df, cbar=0, linewidths=2, square=True, cmap='coolwarm')
plt.savefig('%sradial_heatmap.pdf' % output_dir)
plt.show()