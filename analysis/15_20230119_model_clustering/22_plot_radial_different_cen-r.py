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
data_dir = "%stxt/dataset3/" % master_folder
output_dir = "%sfigures1/" % master_folder

ac = 5000

seps = [0, 10, 20, 30, 40, 50, 60, 70]
sep_feature = 'cen_r'
cmap = mpl.cm.get_cmap('Spectral')
x = np.arange(0, 1, 1/len(seps))
line_color = [cmap(i) for i in x]

data = pd.DataFrame()
for j in seps:
    df = pd.read_csv(("%s%s/%s_5_cluster.txt" % (data_dir, j, ac)), na_values=['.'], sep='\t')
    feature = ['radial_curve_DNAFISH']
    for f in feature:
        df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]
    df['r'] = np.sqrt(df['area_nuclear']/math.pi)
    data = pd.concat([data, df], axis=0)

# radial curve
x = np.arange(0.025, 1, 0.05)

mean_df = []
plt.subplots(figsize=(12, 9))
for i in range(len(seps)):
    data_sample = data[data[sep_feature] == seps[i]].copy().reset_index(drop=True)
    for j in range(len(data_sample)):
        plt.plot(x, data_sample['radial_curve_DNAFISH'][j], alpha=0.01, color=line_color[i])
    mean_curve, ci_lower, ci_higher = dat.mean_list(data_sample['radial_curve_DNAFISH'].tolist())
    mean_df.append(mean_curve)
    plt.plot(x, mean_curve, color=line_color[i], label='%s, n=%s' % (seps[i], len(data_sample)))
plt.xlabel('relative_r')
plt.ylabel('radial_curve')
plt.ylim([0, 10])
plt.legend(loc=2)
plt.savefig('%sac%s_radial_curve.pdf' % (output_dir, ac))
plt.show()

# heatmap
df = pd.DataFrame(mean_df)
df.columns = ['%.3f' % elem for elem in x]
df.index = seps
plt.subplots(figsize=(12, 9))
sns.heatmap(df, cbar=0, linewidths=2, square=True, cmap='coolwarm')
plt.savefig('%sac%s_radial_heatmap.pdf' % (output_dir, ac))
plt.show()