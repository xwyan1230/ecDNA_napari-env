import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
import skimage.io as skio
import shared.display as dis
import shared.dataframe as dat
import shared.image as ima
import numpy as np
import pandas as pd
import shared.math as mat
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import scipy.cluster.hierarchy as shc
import matplotlib.pyplot as plt
import matplotlib.cm as pcm
import napari
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220826_neuroblastoma/"
group = 'ABC'
group_lst = ['A', 'B', 'C']
save_path = '%sv1_img/%s/' % (master_folder, group)
if not os.path.exists(save_path):
    os.makedirs(save_path)
version = 1
cmap = matplotlib.cm.get_cmap('Spectral')

data_mean = pd.read_csv('%s%s_v%s_mean.txt' % (master_folder, group, version), na_values=['.'], sep='\t')
data_mean['area_ratio_per_ecDNA'] = data_mean['total_area_ratio_ecDNA']/data_mean['n_ecDNA']
data_mean['corrected_area_ratio'] = data_mean['total_area_ratio_ecDNA']/data_mean['dis_to_hub_area']
line_colors = []
for i in np.arange(0, 1, 1/len(group_lst)):
    line_colors.append(cmap(i))
line_colors.append((0.30, 0.30, 0.30, 1.0))

# data_mean_A = data_mean[data_mean['group'] == 'A']
# data_mean_B = data_mean[data_mean['group'] == 'B']
# data_mean_C = data_mean[data_mean['group'] == 'C']
A_list = []
B_list = []
C_list = []

for i in range(50):

    data_mean_A = data_mean[(data_mean['group'] == 'A')].copy().reset_index(drop=True)
    data_mean_B = data_mean[(data_mean['group'] == 'B')].sample(n=14).copy().reset_index(drop=True)
    data_mean_C = data_mean[(data_mean['group'] == 'C')].sample(n=14).copy().reset_index(drop=True)

    plt.figure(figsize=(12, 8))
    x = 'dis_to_hub_area'
    y = 'total_area_ratio_ecDNA'
    plt.scatter(data_mean_A[x], data_mean_A[y], c=line_colors[0])
    plt.scatter(data_mean_B[x], data_mean_B[y], c=line_colors[1])
    plt.scatter(data_mean_C[x], data_mean_C[y], c=line_colors[2])

    _, int_fit_r2_sample, int_fit_a_sample = mat.fitting_linear_b0(data_mean['dis_to_hub_area'].tolist(), data_mean['total_area_ratio_ecDNA'].tolist())
    # print(int_fit_a_sample)

    A_above = len(data_mean_A[(data_mean_A['total_area_ratio_ecDNA'] > 0.05) & (data_mean_A['corrected_area_ratio'] > 0.002455)])
    B_above = len(data_mean_B[(data_mean_B['total_area_ratio_ecDNA'] > 0.05) & (data_mean_B['corrected_area_ratio'] > 0.002455)])
    C_above = len(data_mean_C[(data_mean_C['total_area_ratio_ecDNA'] > 0.05) & (data_mean_C['corrected_area_ratio'] > 0.002455)])
    A_list.append(A_above/14.0)
    B_list.append(B_above/14.0)
    C_list.append(C_above/14.0)

    """x = np.arange(0, 50, 5)
    y_sample = int_fit_a_sample * x
    # b = 2 * int_fit_a_sample / (1 - int_fit_a_sample ** 2)
    # y_sample1 = b * x
    plt.plot(x, y_sample, linewidth=2, linestyle='--')
    plt.axhline(0.05, linestyle='--')
    # plt.plot(x, y_sample1, linewidth=2, linestyle='--')
    plt.close()"""

print(A_list)
print(B_list)
print(C_list)
print(np.sum(A_list)/50)
print(np.sum(B_list)/50)
print(np.sum(C_list)/50)
