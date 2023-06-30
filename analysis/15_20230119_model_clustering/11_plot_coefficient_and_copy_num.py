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
import math
import os
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230119_model_clustering/"
data_dir = "%stxt/dataset4/" % master_folder
output_dir = "%sfigures4/" % master_folder

samples = ['0_5', '1_5', '2_5', '5_5', '10_5', '25_5', '50_5', '75_5', '100_5', '200_5', '300_5', '400_5', '500_5', '1000_5',
           '2000_5', '3000_5', '4000_5', '5000_5']
interval = 25
copy_nums = np.arange(0, 201, interval)

g1 = [[] for _ in range(len(samples))]
g20 = [[] for _ in range(len(samples))]
g40 = [[] for _ in range(len(samples))]
above_area = [[] for _ in range(len(samples))]
k_slope = [[] for _ in range(len(samples))]
b_intercept = [[] for _ in range(len(samples))]
zero_percent = [[] for _ in range(len(samples))]
total_area = [[] for _ in range(len(samples))]
dis_to_hub_area = [[] for _ in range(len(samples))]

for j in range(len(samples)):
    data = pd.read_csv(("%s%s_cluster.txt" % (data_dir, samples[j])), na_values=['.'], sep='\t')
    data['r'] = np.sqrt(data['area_nuclear'] / math.pi)
    data['total_area_ecDNA_sqrt'] = np.sqrt(data['total_area_ecDNA'] / math.pi)
    data['total_area_ecDNA_sqrt_normalized'] = data['total_area_ecDNA_sqrt'] / data['r']
    data['dis_to_hub_area_v2_normalized'] = data['dis_to_hub_area_v2'] / data['r']
    feature = ['cum_area_ind_ecDNA_filled', 'percentage_area_curve_ecDNA', 'g']
    for f in feature:
        data[f] = [dat.str_to_float(data[f][i]) for i in range(len(data))]
    data['total_area_ecDNA_sqrt'] = np.sqrt(data['total_area_ecDNA'])
    for i in range(len(copy_nums)-1):
        coefficient = int(samples[j].split('_')[0])
        data_sample = data[(data['copy_num'] > copy_nums[i]) & (data['copy_num'] <= copy_nums[i+1])].copy().reset_index(drop=True)
        dis_to_hub_area[j].append(np.mean(data_sample['dis_to_hub_area_v2_normalized']))
        mean_curve, ci_lower, ci_higher = dat.mean_list(data_sample['g'].tolist())
        print(mean_curve[1])
        g1[j].append(mean_curve[1])
        g20[j].append(mean_curve[20])
        g40[j].append(mean_curve[40])
        mean_curve, ci_lower, ci_higher = dat.mean_list(dat.list_fill_with_last_num(data_sample['percentage_area_curve_ecDNA'].tolist()))
        above_area_temp = len(mean_curve) - 0.5 - sum(mean_curve)
        above_area[j].append(above_area_temp)
        mean_curve, ci_lower, ci_higher = dat.mean_list(dat.list_fill_with_last_num(data_sample['cum_area_ind_ecDNA_filled'].tolist()))
        total_area[j].append(mean_curve[-1])
        x = np.array(data_sample['total_area_ecDNA_sqrt']).reshape((-1, 1))
        y = np.array(data_sample['dis_to_hub_area_v2'])
        model = LinearRegression().fit(x, y)
        print('sample: %s, %s, slope: %s, intercept: %s, r_square: %s' % (coefficient, copy_nums[i+1], model.coef_[0],
                                                                          model.intercept_, model.score(x, y)))
        k_slope[j].append(model.coef_[0])
        b_intercept[j].append(model.intercept_)
        zero_percent[j].append(len(data_sample[data_sample['dis_to_hub_area_v2'] == 0]) * 1.0 / len(data_sample))

df = pd.DataFrame(dis_to_hub_area)
df.columns = copy_nums[1:]
df.index = samples
plt.subplots(figsize=(12, 9))
sns.heatmap(df)
plt.savefig('%s2d_dis_to_hub_area.pdf' % output_dir)
plt.show()

df = pd.DataFrame(g1)
df.columns = copy_nums[1:]
df.index = samples
plt.subplots(figsize=(12, 9))
sns.heatmap(df)
plt.savefig('%s2d_g1.pdf' % output_dir)
plt.show()

df = pd.DataFrame(g20)
df.columns = copy_nums[1:]
df.index = samples
plt.subplots(figsize=(12, 9))
sns.heatmap(df)
plt.savefig('%s2d_g20.pdf' % output_dir)
plt.show()

df = pd.DataFrame(g40)
df.columns = copy_nums[1:]
df.index = samples
plt.subplots(figsize=(12, 9))
sns.heatmap(df)
plt.savefig('%s2d_g40.pdf' % output_dir)
plt.show()

df = pd.DataFrame(above_area)
df.columns = copy_nums[1:]
df.index = samples
plt.subplots(figsize=(12, 9))
sns.heatmap(df)
plt.savefig('%s2d_above_area.pdf' % output_dir)
plt.show()

df = pd.DataFrame(k_slope)
df.columns = copy_nums[1:]
df.index = samples
plt.subplots(figsize=(12, 9))
sns.heatmap(df)
plt.savefig('%s2d_k_slope.pdf' % output_dir)
plt.show()

df = pd.DataFrame(b_intercept)
df.columns = copy_nums[1:]
df.index = samples
plt.subplots(figsize=(12, 9))
sns.heatmap(df)
plt.savefig('%s2d_b_intercept.pdf' % output_dir)
plt.show()

df = pd.DataFrame(zero_percent)
df.columns = copy_nums[1:]
df.index = samples
plt.subplots(figsize=(12, 9))
sns.heatmap(df)
plt.savefig('%s2d_zero_percent.pdf' % output_dir)
plt.show()

df = pd.DataFrame(total_area)
df.columns = copy_nums[1:]
df.index = samples
plt.subplots(figsize=(12, 9))
sns.heatmap(df)
plt.savefig('%s2d_total_area.pdf' % output_dir)
plt.show()
