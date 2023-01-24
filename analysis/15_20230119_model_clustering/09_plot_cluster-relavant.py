import skimage.io as skio
import pandas as pd
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.math as mat
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import seaborn as sns
import shared.dataframe as dat
import os
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230119_model_clustering/"
data_dir = "%stxt/" % master_folder
output_dir = "%sfigures3/" % master_folder

samples = ['0_5', '1_5', '2_5', '5_5', '10_5', '25_5', '50_5', '75_5', '100_5', '200_5', '300_5', '400_5', '500_5', '1000_5',
           '2000_5', '3000_5', '4000_5', '5000_5']

cmap = mpl.cm.get_cmap('Spectral')
x = np.arange(0, 1, 1/len(samples))
line_color = [cmap(i) for i in x]

data = pd.DataFrame()
for i in samples:
    df = pd.read_csv(("%s%s_cluster.txt" % (data_dir, i)), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0)

data['total_area_ecDNA_sqrt'] = np.sqrt(data['total_area_ecDNA'])

# plot zero percentage for all the samples
zero_percentage = []
plt.subplots(figsize=(12, 9))
for i in range(len(samples)):
    coefficient = int(samples[i].split('_')[0])
    data_sample = data[data['coefficient'] == coefficient].copy().reset_index(drop=True)
    zero_percentage.append(len(data_sample[data_sample['dis_to_hub_area_v2'] == 0])*1.0/len(data_sample))
plt.bar(samples, zero_percentage, color=line_color)
plt.xlabel('sample coefficient')
plt.ylabel('percentage of zero dis_to_hub_area_v2')
plt.savefig('%szero_percentage.pdf' % output_dir)
plt.show()

# plot total area for all the samples
sns.set_palette(sns.color_palette(line_color))
plt.subplots(figsize=(12, 9))
sns.violinplot(data=data, x='coefficient', y='total_area_ecDNA', hue_order=samples)
plt.savefig('%stotal_area_ecDNA.pdf' % output_dir)
plt.show()

# plot linear fitting of non-zero measurements a b for all the samples
slopes = []
intercepts = []
r_sq = []
for i in range(len(samples)):
    coefficient = int(samples[i].split('_')[0])
    data_sample = data[data['coefficient'] == coefficient].copy().reset_index(drop=True)
    x = np.array(data_sample['total_area_ecDNA_sqrt']).reshape((-1, 1))
    y = np.array(data_sample['dis_to_hub_area_v2'])
    model = LinearRegression().fit(x, y)
    print('coefficient: %s, slope: %s, intercept: %s, r_square: %s' % (coefficient, model.coef_[0], model.intercept_,
                                                                       model.score(x, y)))
    slopes.append(model.coef_[0])
    intercepts.append(model.intercept_)
    r_sq.append(model.score(x, y))

plt.subplots(figsize=(12, 9))
plt.bar(samples, slopes, color=line_color)
plt.savefig('%slinearregression_slope.pdf' % output_dir)
plt.xlabel('sample coefficient')
plt.ylabel('linear regression slope')
plt.show()

plt.subplots(figsize=(12, 9))
plt.bar(samples, intercepts, color=line_color)
plt.savefig('%slinearregression_intercept.pdf' % output_dir)
plt.xlabel('sample coefficient')
plt.ylabel('linear regression intercept')
plt.show()

plt.subplots(figsize=(12, 9))
plt.bar(samples, r_sq, color=line_color)
plt.savefig('%slinearregression_rsq.pdf' % output_dir)
plt.xlabel('sample coefficient')
plt.ylabel('linear regression r square')
plt.show()