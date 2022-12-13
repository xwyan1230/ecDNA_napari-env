import matplotlib.pyplot as plt
import skimage.io as skio
import napari
import pandas as pd
import shared.dataframe as dat
import seaborn as sns
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221208_analysis_20221130_EdU_ON_metaphase/figures/"
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221208_analysis_20221130_EdU_ON_metaphase/figures/"
output_dir = data_dir
sample = 'point5uM_EdU'

# load and format data
df = pd.read_csv(("%s%s.txt" % (master_folder, sample)), na_values=['.'], sep='\t')
feature = ['EdU_ind_mean_int', 'EdU_ind_area', 'MYC_ind_mean_int']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]
df['EdU_centroid'] = [dat.str_to_list_of_float(df['EdU_centroid'][i], 2) for i in range(len(df))]
df['EdU_mean'] = [np.mean(df['EdU_ind_mean_int'][i]) for i in range(len(df))]
df['MYC_mean'] = [np.mean(df['MYC_ind_mean_int'][i]) for i in range(len(df))]

df = df.sort_values(by='EdU_chr_mean_int').copy().reset_index(drop=True)

index_temp = []
fov_temp = []
nuclear_temp = []
EdU_mean_int_temp = []
MYC_mean_int_temp = []
EdU_chr_mean_int_temp = []
for i in range(len(df)):
    index_temp = index_temp + [df.index[i]] * df['n'][i]
    fov_temp = fov_temp + [df['FOV'][i]] * df['n'][i]
    nuclear_temp = nuclear_temp + [df['nuclear'][i]] * df['n'][i]
    EdU_mean_int_temp = EdU_mean_int_temp + df['EdU_ind_mean_int'][i]
    MYC_mean_int_temp = MYC_mean_int_temp + df['MYC_ind_mean_int'][i]
    EdU_chr_mean_int_temp = EdU_chr_mean_int_temp + [df['EdU_chr_mean_int'][i]] * df['n'][i]

EdU = pd.DataFrame({'index': index_temp, 'FOV': fov_temp, 'nuclear': nuclear_temp, 'EdU_mean_int': EdU_mean_int_temp,
                    'MYC_mean_int': MYC_mean_int_temp, 'EdU_chr_mean_int': EdU_chr_mean_int_temp})
EdU['EdU_normalized_mean_int'] = EdU['EdU_mean_int']/EdU['MYC_mean_int']

fig, ax = plt.subplots(figsize=(9,7))
sp = sns.scatterplot(data=EdU, x='MYC_mean_int', y='EdU_mean_int', c=EdU['EdU_chr_mean_int'], s=7, alpha=1, label='cell=%s, n=%s' % (len(df), len(EdU)))
plt.legend()
plt.savefig('%s%s_scatter_EdU_vs_MYC.tiff' % (output_dir, sample))
plt.close()

for i in range(len(df)):
    fov = df['FOV'][i]
    nuclear = df['nuclear'][i]
    sns.scatterplot(data=EdU, x='MYC_mean_int', y='EdU_mean_int', color='black', s=5, alpha=0.1)
    sns.scatterplot(data=EdU[EdU['index'] == i], x='MYC_mean_int', y='EdU_mean_int', color='r', s=5, alpha=1)
    plt.savefig('%s%s_%s_%s_scatter_EdU_vs_MYC.tiff' % (output_dir, sample, fov, nuclear))
    plt.close()

