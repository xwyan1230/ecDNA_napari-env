import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import shared.dataframe as dat
import seaborn as sns
from shared.sinaplot import sinaplot
from scipy.stats import ks_2samp
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230602_analysis_PC3DMandHSR_1uM_6hr_oldMYC/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_nuclear_convex_dilation = -3

# samples
exp1 = 'PC3HSR_6hr'
exp2 = 'PC3DM_6hr'
ctrl1 = pd.read_csv('%s%s_ctrl.txt' % (data_dir, exp1), na_values=['.'], sep='\t')
ctrl2 = pd.read_csv('%s%s_ctrl.txt' % (data_dir, exp2), na_values=['.'], sep='\t')
ctrl1['well'] = [ctrl1['sample'][i].split('_')[0] for i in range(len(ctrl1))]
ctrl2['well'] = [ctrl2['sample'][i].split('_')[0] for i in range(len(ctrl2))]
feature1 = ['total_int_MYC_mean_log2FC', 'mean_int_MYC_mean_log2FC']
feature2 = ['total_int_MYC_mean_log2FC', 'mean_int_MYC_mean_log2FC']

# load data
df1 = pd.read_csv('%s%s_n%s_sum.txt' % (data_dir1, exp1, n_nuclear_convex_dilation), na_values=['.'], sep='\t')
df2 = pd.read_csv('%s%s_n%s_sum.txt' % (data_dir1, exp2, n_nuclear_convex_dilation), na_values=['.'], sep='\t')

df_WT1 = pd.DataFrame()
df_WT2 = pd.DataFrame()
for i in range(len(ctrl1)):
    df_temp1 = df1[(df1['plate'] == ctrl1['plate'][i]) & (df1['well'] == ctrl1['well'][i])].copy().reset_index(drop=True)
    df_temp2 = df2[(df2['plate'] == ctrl2['plate'][i]) & (df2['well'] == ctrl2['well'][i])].copy().reset_index(drop=True)
    df_WT1 = pd.concat([df_WT1, df_temp1], axis=0).reset_index(drop=True)
    df_WT2 = pd.concat([df_WT2, df_temp2], axis=0).reset_index(drop=True)

# scatter plot
for f in range(len(feature1)):
    model = LinearRegression().fit(np.array(df1[feature1[f]]).reshape((-1, 1)), np.array(df2[feature2[f]]))
    b = model.intercept_
    a = model.coef_[0]
    plt.subplots(figsize=(8, 8))
    plt.scatter(df1[feature1[f]], df2[feature2[f]], color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df1)):
        plt.text(x=df1[feature1[f]][i]+0.005, y=df2[feature2[f]][i]+0.005, s=df1['p-w'][i], size=5)
    plt.scatter(df_WT1[feature1[f]], df_WT2[feature2[f]], color=(0.85, 0.35, 0.25), s=20)
    limit = max(abs(df1[feature1[f]].max()), abs(df1[feature1[f]].min()), abs(df2[feature2[f]].max()),
                abs(df2[feature2[f]].min()))
    x = np.arange(-limit, limit, 0.01)
    y = a*x+b
    plt.plot(x, y, linestyle='--', color=(0.85, 0.35, 0.25))
    plt.xlabel(exp1)
    plt.ylabel(exp2)
    """limit1 = max(abs(df1[feature1[f]].max()), abs(df1[feature1[f]].min()))
    plt.xlim([-(limit1+0.01), (limit1+0.01)])
    limit2 = max(abs(df2[feature2[f]].max()), abs(df2[feature2[f]].min()))
    plt.ylim([-(limit2 + 0.01), (limit2 + 0.01)])"""
    plt.xlim([-(limit + 0.01), (limit + 0.01)])
    plt.ylim([-(limit + 0.01), (limit + 0.01)])
    plt.savefig('%s/%s_vs_%s_%s_pw_linear.pdf' % (output_dir, exp1, exp2, feature1[f]))
    plt.close()

    plt.subplots(figsize=(8, 8))
    plt.scatter(df1[feature1[f]], df2[feature2[f]], color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df1)):
        plt.text(x=df1[feature1[f]][i] + 0.005, y=df2[feature2[f]][i] + 0.005, s=df1['compound'][i], size=5)
    plt.scatter(df_WT1[feature1[f]], df_WT2[feature2[f]], color=(0.85, 0.35, 0.25), s=20)
    limit = max(abs(df1[feature1[f]].max()), abs(df1[feature1[f]].min()), abs(df2[feature2[f]].max()),
                abs(df2[feature2[f]].min()))
    x = np.arange(-limit, limit, 0.01)
    y = a * x + b
    plt.plot(x, y, linestyle='--', color=(0.85, 0.35, 0.25))
    plt.xlabel(exp1)
    plt.ylabel(exp2)
    """limit1 = max(abs(df1[feature1[f]].max()), abs(df1[feature1[f]].min()))
    plt.xlim([-(limit1+0.01), (limit1+0.01)])
    limit2 = max(abs(df2[feature2[f]].max()), abs(df2[feature2[f]].min()))
    plt.ylim([-(limit2 + 0.01), (limit2 + 0.01)])"""
    plt.xlim([-(limit + 0.01), (limit + 0.01)])
    plt.ylim([-(limit + 0.01), (limit + 0.01)])
    plt.savefig('%s/%s_vs_%s_%s_compound_linear.pdf' % (output_dir, exp1, exp2, feature1[f]))
    plt.close()

    plt.subplots(figsize=(8, 8))
    plt.scatter(df1[feature1[f]], df2[feature2[f]], color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df1)):
        plt.text(x=df1[feature1[f]][i] + 0.005, y=df2[feature2[f]][i] + 0.005, s=df1['target'][i], size=5)
    plt.scatter(df_WT1[feature1[f]], df_WT2[feature2[f]], color=(0.85, 0.35, 0.25), s=20)
    limit = max(abs(df1[feature1[f]].max()), abs(df1[feature1[f]].min()), abs(df2[feature2[f]].max()),
                abs(df2[feature2[f]].min()))
    x = np.arange(-limit, limit, 0.01)
    y = a * x + b
    plt.plot(x, y, linestyle='--', color=(0.85, 0.35, 0.25))
    plt.xlabel(exp1)
    plt.ylabel(exp2)
    """limit1 = max(abs(df1[feature1[f]].max()), abs(df1[feature1[f]].min()))
    plt.xlim([-(limit1+0.01), (limit1+0.01)])
    limit2 = max(abs(df2[feature2[f]].max()), abs(df2[feature2[f]].min()))
    plt.ylim([-(limit2 + 0.01), (limit2 + 0.01)])"""
    plt.xlim([-(limit + 0.01), (limit + 0.01)])
    plt.ylim([-(limit + 0.01), (limit + 0.01)])
    plt.savefig('%s/%s_vs_%s_%s_target_linear.pdf' % (output_dir, exp1, exp2, feature1[f]))
    plt.close()