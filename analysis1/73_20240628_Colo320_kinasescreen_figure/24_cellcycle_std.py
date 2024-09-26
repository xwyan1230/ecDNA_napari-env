import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import scipy.cluster.hierarchy as shc
from pandas.plotting import parallel_coordinates
import matplotlib.pyplot as plt
import matplotlib.cm as pcm
import matplotlib
import seaborn as sns
import shared.dataframe as dat
import shared.display as dis
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']

df = pd.read_csv('%s/cell-cycle_sampling.txt' % (data_dir), na_values=['.'], sep='\t')
sample_lst = np.arange(1, 3000, 1)

data = pd.DataFrame(columns=['sample_n', 'G1_mean', 'G1_std', 'S_mean', 'S_std', 'G2M_mean', 'G2M_std'])
for s in sample_lst:
    data_temp = df[df['sample_n']==s].copy().reset_index(drop=True)
    data.loc[len(data.index)] = [s, np.mean(data_temp['per_G1']), np.std(data_temp['per_G1']),
                                 np.mean(data_temp['per_S']), np.std(data_temp['per_S']),
                                 np.mean(data_temp['per_G2M']), np.std(data_temp['per_G2M'])]

data.to_csv('%s/cell-cycle_std.txt' % (data_dir), index=False, sep='\t')


def model_func(x, k, b):
    return np.exp(k*np.log(x)+b)


def model_func1(x, k):
    return k


features = ['G1_mean', 'G1_std', 'S_mean', 'S_std', 'G2M_mean', 'G2M_std']
for feature in features:
    print(feature)
    plt.subplots(figsize=(9, 7))
    sns.scatterplot(data=data, x='sample_n', y=feature, s=5, color=rainboo_colors[8])
    if feature[-1] == 'd':
        popt, _ = curve_fit(model_func, data['sample_n'].tolist(), data[feature].tolist())
        k, b = popt
        print(k)
        print(b)
        fit_y = model_func(np.array(sample_lst), k, b)
        plt.plot(sample_lst, fit_y, '--', color=rainboo_colors[1])
    else:
        popt, _ = curve_fit(model_func1, data['sample_n'].tolist(), data[feature].tolist())
        k = popt[0]
        print(k)
        plt.plot(sample_lst, [k] * len(sample_lst), '--', color=rainboo_colors[1])
        # plt.axhline(y=k, linestyle='--', color=rainboo_colors[1])
    # plt.ylim([0, 1.05])
    plt.savefig('%s/mix_theoretical_%s.pdf' % (output_dir, feature))
    plt.show()
    """if feature[-1] == 'd':
        plt.subplots(figsize=(9, 7))
        data['temp_x'] = np.log(data['sample_n'])
        data['temp_y'] = np.log(data[feature])
        sns.scatterplot(data=data, x='temp_x', y='temp_y', s=5)
        plt.show()"""



