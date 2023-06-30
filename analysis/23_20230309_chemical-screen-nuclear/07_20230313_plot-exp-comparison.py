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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230309_analysis_chemical-screen-nuclear/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_nuclear_convex_dilation = -3

# samples
exp1 = 'HSR_2hr'
exp2 = 'DM_2hr'
ctrls1 = ['HSR_2hr*C3', 'HSR_2hr*C10', 'HSR_2hr*D6', 'HSR_2hr*F3', 'HSR_2hr*F10']
ctrls2 = ['DM_2hr*C3', 'DM_2hr*C10', 'DM_2hr*D6', 'DM_2hr*F3', 'DM_2hr*F10']
features = ['mean_int_MYC_mean']

# load data
df1 = pd.read_csv('%s%s_n%s_sum.txt' % (data_dir, exp1, n_nuclear_convex_dilation), na_values=['.'], sep='\t')
df2 = pd.read_csv('%s%s_n%s_sum.txt' % (data_dir, exp2, n_nuclear_convex_dilation), na_values=['.'], sep='\t')

df_WT1 = pd.DataFrame()
df_WT2 = pd.DataFrame()
for i in range(len(ctrls1)):
    df_temp1 = df1[(df1['plate'] == ctrls1[i].split('*')[0]) & (df1['well'] == ctrls1[i].split('*')[1])].copy().reset_index(drop=True)
    df_temp2 = df2[(df2['plate'] == ctrls2[i].split('*')[0]) & (df2['well'] == ctrls2[i].split('*')[1])].copy().reset_index(drop=True)
    df_WT1 = pd.concat([df_WT1, df_temp1], axis=0).reset_index(drop=True)
    df_WT2 = pd.concat([df_WT2, df_temp2], axis=0).reset_index(drop=True)

# scatter plot
for feature in features:
    plt.subplots(figsize=(9, 6))
    plt.scatter(df1[feature], df2[feature], color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df1)):
        plt.text(x=df1[feature][i]+0.5, y=df2[feature][i]+0.5, s=df1['p-w'][i], size=5)
    plt.scatter(df_WT1[feature], df_WT2[feature], color=(0.85, 0.35, 0.25), s=20)
    plt.xlabel(exp1)
    plt.ylabel(exp2)
    plt.xlim([0, 33000])
    plt.ylim([0, 33000])
    plt.savefig('%s/%s_vs_%s_%s_pw.pdf' % (output_dir, exp1, exp2, feature))
    plt.close()

    plt.subplots(figsize=(9, 6))
    plt.scatter(df1[feature], df2[feature], color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df1)):
        plt.text(x=df1[feature][i] + 0.5, y=df2[feature][i] + 0.5, s=df1['compound'][i], size=5)
    plt.scatter(df_WT1[feature], df_WT2[feature], color=(0.85, 0.35, 0.25), s=20)
    plt.xlabel(exp1)
    plt.ylabel(exp2)
    plt.xlim([0, 33000])
    plt.ylim([0, 33000])
    plt.savefig('%s/%s_vs_%s_%s_compound.pdf' % (output_dir, exp1, exp2, feature))
    plt.close()

    plt.subplots(figsize=(9, 6))
    plt.scatter(df1[feature], df2[feature], color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df1)):
        plt.text(x=df1[feature][i] + 0.5, y=df2[feature][i] + 0.5, s=df1['target'][i], size=5)
    plt.scatter(df_WT1[feature], df_WT2[feature], color=(0.85, 0.35, 0.25), s=20)
    plt.xlabel(exp1)
    plt.ylabel(exp2)
    plt.xlim([0, 33000])
    plt.ylim([0, 33000])
    plt.savefig('%s/%s_vs_%s_%s_target.pdf' % (output_dir, exp1, exp2, feature))
    plt.close()