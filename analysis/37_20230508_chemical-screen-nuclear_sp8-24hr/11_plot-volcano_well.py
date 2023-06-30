import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
from shared.sinaplot import sinaplot
from scipy.stats import ks_2samp
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230508_analysis_chemical-screen-nuclear_24hr/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_nuclear_convex_dilation = -3

# samples
exp = 'HSR_24hr'
ctrl = pd.read_csv('%s%s_ctrl.txt' % (data_dir, exp), na_values=['.'], sep='\t')
features = ['area_nuclear', 'mean_int_nuclear', 'mean_int_MYC', 'mean_int_nuclear_cal', 'mean_int_MYC_cal',
            'total_int_nuclear', 'total_int_MYC', 'total_int_nuclear_cal', 'total_int_MYC_cal']
ctrl['well'] = [ctrl['sample'][i].split('_')[0] for i in range(len(ctrl))]

# load data
df = pd.read_csv('%s%s_n%s_sum.txt' % (data_dir1, exp, n_nuclear_convex_dilation), na_values=['.'], sep='\t')

df_WT = pd.DataFrame()
for i in range(len(ctrl)):
    df_temp = df[(df['plate'] == ctrl['plate'][i]) & (df['well'] == ctrl['well'][i])].copy().reset_index(drop=True)
    df_WT = pd.concat([df_WT, df_temp], axis=0).reset_index(drop=True)

# scatter plot
for feature in features:
    plt.subplots(figsize=(9, 6))
    sns.scatterplot(data=df, y='%s_p' % feature, x='%s_mean_log2FC' % feature, color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df)):
        plt.text(x=df['%s_mean_log2FC' % feature][i]+0.001, y=df['%s_p' % feature][i]+0.1, s=df['p-w'][i], size=5)
    sns.scatterplot(data=df_WT, y='%s_p' % feature, x='%s_mean_log2FC' % feature, color=(0.85, 0.35, 0.25), s=20)
    plt.axhline(y=-np.log(0.01), linestyle='--', color='red')
    limit = max(abs(df['%s_mean_log2FC' % feature].max()), abs(df['%s_mean_log2FC' % feature].min()))
    plt.xlim([-(limit+0.01), limit+0.01])
    plt.savefig('%s/%s_%s_volcano_pw.pdf' % (output_dir, exp, feature))
    plt.close()

    plt.subplots(figsize=(9, 6))
    sns.scatterplot(data=df, y='%s_p' % feature, x='%s_mean_log2FC' % feature, color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df)):
        plt.text(x=df['%s_mean_log2FC' % feature][i] + 0.001, y=df['%s_p' % feature][i] + 0.1, s=df['compound'][i], size=5)
    sns.scatterplot(data=df_WT, y='%s_p' % feature, x='%s_mean_log2FC' % feature, color=(0.85, 0.35, 0.25), s=20)
    plt.axhline(y=-np.log(0.01), linestyle='--', color='red')
    limit = max(abs(df['%s_mean_log2FC' % feature].max()), abs(df['%s_mean_log2FC' % feature].min()))
    plt.xlim([-(limit + 0.01), limit + 0.01])
    plt.savefig('%s/%s_%s_volcano_compound.pdf' % (output_dir, exp, feature))
    plt.close()

    plt.subplots(figsize=(9, 6))
    sns.scatterplot(data=df, y='%s_p' % feature, x='%s_mean_log2FC' % feature, color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df)):
        plt.text(x=df['%s_mean_log2FC' % feature][i] + 0.001, y=df['%s_p' % feature][i] + 0.1, s=df['target'][i], size=5)
    sns.scatterplot(data=df_WT, y='%s_p' % feature, x='%s_mean_log2FC' % feature, color=(0.85, 0.35, 0.25), s=20)
    plt.axhline(y=-np.log(0.01), linestyle='--', color='red')
    limit = max(abs(df['%s_mean_log2FC' % feature].max()), abs(df['%s_mean_log2FC' % feature].min()))
    plt.xlim([-(limit + 0.01), limit + 0.01])
    plt.savefig('%s/%s_%s_volcano_target.pdf' % (output_dir, exp, feature))
    plt.close()