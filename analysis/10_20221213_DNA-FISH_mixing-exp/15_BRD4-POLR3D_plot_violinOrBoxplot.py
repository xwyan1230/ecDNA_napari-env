import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221213_analysis_DNA-FISH_mixing_exp/20221121_H2B-series_POLR3D-BRD4-BRD1-DAPK2/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'H2B+BRD1'
pos_threshold = 30000
neg_threshold = 20000
pos = 'ctrl'
neg = 'BRD1 KO'
hue_order = [pos, neg]

df = pd.read_csv(("%s%s_ecDNA.txt" % (data_dir, sample)), na_values=['.'], sep='\t')
df['total_int_DNAFISH'] = df['mean_int_DNAFISH'] * df['area_nuclear']
df['total_int_ecDNA'] = df['mean_int_ecDNA'] * df['total_area_ecDNA']

hist_colors = [(0.90, 0.90, 0.90), (0.95, 0.50, 0.50), (0.50, 0.90, 0.90)]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25), (0.30, 0.70, 0.70)]

sample_lst = []
for i in range(len(df)):
    if df['mCherry_mean'].tolist()[i] < neg_threshold:
        sample_lst.append(neg)
    elif df['mCherry_mean'].tolist()[i] > pos_threshold:
        sample_lst.append(pos)
    else:
        sample_lst.append('NA')
df['sample'] = sample_lst

df_sort = df[df['sample'].isin([neg, pos])].copy().reset_index(drop=True)

# plot multiple feature in one graph
"""
feature = ['mean_int_nuclear']
feature_lst = []
sample_lst = []
value_lst = []
data = pd.DataFrame()
for i in range(len(feature)):
    feature_lst = feature_lst + [feature[i]] * len(df_sort)
    sample_lst = sample_lst + df_sort['sample'].tolist()
    value_lst = value_lst + df_sort[feature[i]].tolist()
data['feature'] = feature_lst
data['sample'] = sample_lst
data['value'] = value_lst

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(3*len(feature), 9))
sns.barplot(data=data, x='feature', y='value', estimator=np.mean, hue='sample')
plt.savefig('%s/mean_int_nuclear_%s.pdf' % (output_dir, sample))
plt.show()"""

# generate normalize_df
"""feature = ['mean_int_nuclear', 'area_nuclear', 'n_ecDNA', 'total_int_DNAFISH', 'total_int_ecDNA', 'total_area_ecDNA', 'total_area_ratio_ecDNA']

normalize_df = pd.DataFrame()
normalize_df['sample'] = [neg, pos]

for i in range(len(feature)):
    data = pd.DataFrame()
    data['feature'] = [feature[i]] * len(df_sort)
    data['sample'] = df_sort['sample'].tolist()
    data['value'] = df_sort[feature[i]].tolist()

    normalize_df[feature[i]] = [np.mean(data[data['sample'] == neg]['value']), np.mean(data[data['sample'] == pos]['value'])]

normalize_df.to_csv('%snormalize_df.txt' % output_dir, index=False, sep='\t')"""

# plot feature in absolute value
"""feature = ['mean_int_nuclear', 'area_nuclear', 'n_ecDNA', 'total_int_DNAFISH', 'total_int_ecDNA', 'total_area_ecDNA', 'total_area_ratio_ecDNA']

for i in range(len(feature)):
    data = pd.DataFrame()
    data['feature'] = [feature[i]] * len(df_sort)
    data['sample'] = df_sort['sample'].tolist()
    data['value'] = df_sort[feature[i]].tolist()

    sns.set_palette(sns.color_palette(line_colors))
    fig, ax = plt.subplots(figsize=(3, 9))
    fig.subplots_adjust(left=0.2)
    ax = sns.barplot(data=data, x='feature', y='value', estimator=np.mean, hue='sample', hue_order=hue_order)
    # plt.savefig('%s/barplot_%s_%s.pdf' % (output_dir, feature[i], sample))
    plt.close()
"""
# plot feature in normalized value
feature = ['mean_int_nuclear', 'area_nuclear', 'n_ecDNA', 'total_int_DNAFISH', 'total_int_ecDNA', 'total_area_ecDNA', 'total_area_ratio_ecDNA']
normalize_df = pd.read_csv(("%snormalize_df.txt" % data_dir), na_values=['.'], sep='\t')

for i in range(len(feature)):
    normalize_ratio = normalize_df[normalize_df['sample'] == 'DM H2B-mCherry'][feature[i]].tolist()[0]/normalize_df[normalize_df['sample'] == 'DM'][feature[i]].tolist()[0]
    data = pd.DataFrame()
    data['feature'] = [feature[i]] * len(df_sort)
    data['sample'] = df_sort['sample'].tolist()
    data['value'] = df_sort[feature[i]].tolist()

    data['normalized_value'] = [data['value'][j]/normalize_ratio if data['sample'][j] == pos else data['value'][j] for j in range(len(data))]

    sns.set_palette(sns.color_palette(line_colors))
    fig, ax = plt.subplots(figsize=(3, 9))
    fig.subplots_adjust(left=0.2)
    ax = sns.barplot(data=data, x='feature', y='normalized_value', estimator=np.mean, hue='sample', hue_order=hue_order)
    plt.savefig('%s/normalized_barplot_%s_%s.pdf' % (output_dir, feature[i], sample))
    plt.close()


