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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230316_analysis_chemical-screen-nuclear_rep2_Keyence/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

# samples
exp = 'HSR_2hr'
ctrls = ['HSR_2hr*C3', 'HSR_2hr*C10', 'HSR_2hr*D6', 'HSR_2hr*F3', 'HSR_2hr*F10']
features = ['area_nuclear', 'mean_int_nuclear', 'mean_int_MYC']

# load data
df_full = pd.read_csv('%s%s_sum.txt' % (data_dir, exp), na_values=['.'], sep='\t')
df = df_full[df_full['n_total'] >= 100].copy().reset_index(drop=True)

df_WT = pd.DataFrame()
for i in range(len(ctrls)):
    df_temp = df[(df['plate'] == ctrls[i].split('*')[0]) & (df['well'] == ctrls[i].split('*')[1])].copy().reset_index(drop=True)
    df_WT = pd.concat([df_WT, df_temp], axis=0).reset_index(drop=True)

# scatter plot
for feature in features:
    plt.subplots(figsize=(9, 6))
    sns.scatterplot(data=df, y='%s_p' % feature, x='%s_mean' % feature, color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df)):
        plt.text(x=df['%s_mean' % feature][i]+0.1, y=df['%s_p' % feature][i]+0.1, s=df['p-w'][i], size=5)
    sns.scatterplot(data=df_WT, y='%s_p' % feature, x='%s_mean' % feature, color=(0.85, 0.35, 0.25), s=20)
    plt.savefig('%s/%s_%s_volcano_pw.pdf' % (output_dir, exp, feature))
    plt.close()

    plt.subplots(figsize=(9, 6))
    sns.scatterplot(data=df, y='%s_p' % feature, x='%s_mean' % feature, color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df)):
        plt.text(x=df['%s_mean' % feature][i] + 0.1, y=df['%s_p' % feature][i] + 0.1, s=df['compound'][i], size=5)
    sns.scatterplot(data=df_WT, y='%s_p' % feature, x='%s_mean' % feature, color=(0.85, 0.35, 0.25), s=20)
    plt.savefig('%s/%s_%s_volcano_compound.pdf' % (output_dir, exp, feature))
    plt.close()

    plt.subplots(figsize=(9, 6))
    sns.scatterplot(data=df, y='%s_p' % feature, x='%s_mean' % feature, color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df)):
        plt.text(x=df['%s_mean' % feature][i] + 0.1, y=df['%s_p' % feature][i] + 0.1, s=df['target'][i], size=5)
    sns.scatterplot(data=df_WT, y='%s_p' % feature, x='%s_mean' % feature, color=(0.85, 0.35, 0.25), s=20)
    plt.savefig('%s/%s_%s_volcano_target.pdf' % (output_dir, exp, feature))
    plt.close()