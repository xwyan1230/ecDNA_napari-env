import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
from shared.sinaplot import sinaplot
from scipy.stats import ks_2samp
import os
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

# samples
ctrl = pd.read_csv('%sCtrl.txt' % data_dir, na_values=['.'], sep='\t')
features = ['mean_area_nuclear',  'mean_r10', 'mean_r16', 'mean_n_ecDNA']
features1 = ['mean_r10_point2', 'mean_r16_point2', 'mean_averageD_point2', 'std_averageD_point2', 'mean_n_ecDNA_point2', 'per_n_ecDNA10_point2']
features2 = ['mean_total_area_ecDNA', 'mean_total_area_ratio_ecDNA', 'mean_mean_int_DNAFISH']

ex_list = ['D9']
# load data
df_mean = pd.read_csv('%ssummary_mean.txt' % data_dir, na_values=['.'], sep='\t')
df_p = pd.read_csv('%ssummary_p.txt' % data_dir, na_values=['.'], sep='\t')
df_mean_ctrl = df_mean[df_mean['sample'].isin(ctrl['sample'].tolist())].copy().reset_index(drop=True)
df_p_ctrl = df_p[df_p['sample'].isin(ctrl['sample'].tolist())].copy().reset_index(drop=True)
df_mean_ex = df_mean[~df_mean['sample'].isin(ex_list)].copy().reset_index(drop=True)
df_p_ex = df_p[~df_p['sample'].isin(ex_list)].copy().reset_index(drop=True)
df_mean_ctrl_ex = df_mean_ctrl[~df_mean_ctrl['sample'].isin(ex_list)].copy().reset_index(drop=True)
df_p_ctrl_ex = df_p_ctrl[~df_p_ctrl['sample'].isin(ex_list)].copy().reset_index(drop=True)

if not os.path.exists("%s/volcano/" % output_dir):
    os.makedirs("%s/volcano/" % output_dir)
# scatter plot
for feature in features:
    plt.subplots(figsize=(9, 6))
    plt.scatter(df_mean[feature], df_p[feature], color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df_mean)):
        plt.text(x=df_mean[feature][i]+0.001, y=df_p[feature][i]+0.01, s='%s(%s)' % (df_mean['sample'][i], df_mean['n'][i]), size=5)
    plt.scatter(df_mean_ctrl[feature], df_p_ctrl[feature], color=(0.85, 0.35, 0.25), s=20)
    # plt.axhline(y=-np.log(0.01), linestyle='--', color='red')
    limit = max(abs(df_mean[feature].max()), abs(df_mean[feature].min()))
    plt.xlim([-(limit+0.01), limit+0.01])
    plt.xlabel(feature)
    plt.ylabel('-lnp')
    plt.savefig('%s/volcano/%s_volcano_pw.pdf' % (output_dir, feature))
    plt.close()

for feature in features:
    plt.subplots(figsize=(9, 6))
    plt.scatter(df_mean[feature], df_p[feature], color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df_mean)):
        plt.text(x=df_mean[feature][i]+0.001, y=df_p[feature][i]+0.01, s='%s(%s)' % (df_mean['gene'][i], df_mean['n'][i]), size=5)
    plt.scatter(df_mean_ctrl[feature], df_p_ctrl[feature], color=(0.85, 0.35, 0.25), s=20)
    # plt.axhline(y=-np.log(0.01), linestyle='--', color='red')
    limit = max(abs(df_mean[feature].max()), abs(df_mean[feature].min()))
    plt.xlim([-(limit+0.01), limit+0.01])
    plt.xlabel(feature)
    plt.ylabel('-lnp')
    plt.savefig('%s/volcano/%s_volcano_gene.pdf' % (output_dir, feature))
    plt.close()

for feature in features1:
    plt.subplots(figsize=(9, 6))
    plt.scatter(df_mean[feature], df_p[feature], color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df_mean)):
        plt.text(x=df_mean[feature][i]+0.001, y=df_p[feature][i]+0.01, s='%s(%s)' % (df_mean['sample'][i], df_mean['n_point2'][i]), size=5)
    plt.scatter(df_mean_ctrl[feature], df_p_ctrl[feature], color=(0.85, 0.35, 0.25), s=20)
    # plt.axhline(y=-np.log(0.01), linestyle='--', color='red')
    limit = max(abs(df_mean[feature].max()), abs(df_mean[feature].min()))
    plt.xlim([-(limit+0.01), limit+0.01])
    plt.xlabel(feature)
    plt.ylabel('-lnp')
    plt.savefig('%s/volcano/%s_volcano_pw.pdf' % (output_dir, feature))
    plt.close()

for feature in features1:
    plt.subplots(figsize=(9, 6))
    plt.scatter(df_mean[feature], df_p[feature], color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df_mean)):
        plt.text(x=df_mean[feature][i]+0.001, y=df_p[feature][i]+0.01, s='%s(%s)' % (df_mean['gene'][i], df_mean['n_point2'][i]), size=5)
    plt.scatter(df_mean_ctrl[feature], df_p_ctrl[feature], color=(0.85, 0.35, 0.25), s=20)
    # plt.axhline(y=-np.log(0.01), linestyle='--', color='red')
    limit = max(abs(df_mean[feature].max()), abs(df_mean[feature].min()))
    plt.xlim([-(limit+0.01), limit+0.01])
    plt.xlabel(feature)
    plt.ylabel('-lnp')
    plt.savefig('%s/volcano/%s_volcano_gene.pdf' % (output_dir, feature))
    plt.close()

feature = 'per_n_ecDNA10_point2'
plt.subplots(figsize=(9, 6))
plt.scatter(df_mean[feature], df_p[feature], color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df_mean)):
    plt.text(x=df_mean[feature][i]+0.001, y=df_p[feature][i]+0.01, s=df_mean['sample'][i], size=5)
plt.scatter(df_mean_ctrl[feature], df_p_ctrl[feature], color=(0.85, 0.35, 0.25), s=20)
plt.xlim([-3, 3])
plt.xlabel(feature)
plt.ylabel('-lnp')
plt.savefig('%s/volcano/%s_volcano_pw_rescale.pdf' % (output_dir, feature))
plt.close()

plt.subplots(figsize=(9, 6))
plt.scatter(df_mean[feature], df_p[feature], color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df_mean)):
    plt.text(x=df_mean[feature][i]+0.001, y=df_p[feature][i]+0.01, s=df_mean['gene'][i], size=5)
plt.scatter(df_mean_ctrl[feature], df_p_ctrl[feature], color=(0.85, 0.35, 0.25), s=20)
plt.xlim([-3, 3])
plt.xlabel(feature)
plt.ylabel('-lnp')
plt.savefig('%s/volcano/%s_volcano_gene_rescale.pdf' % (output_dir, feature))
plt.close()

feature = 'mean_r10'
plt.subplots(figsize=(9, 6))
plt.scatter(df_mean[feature], df_p[feature], color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df_mean)):
    plt.text(x=df_mean[feature][i]+0.001, y=df_p[feature][i]+0.01, s=df_mean['sample'][i], size=5)
plt.scatter(df_mean_ctrl[feature], df_p_ctrl[feature], color=(0.85, 0.35, 0.25), s=20)
plt.xlim([-0.1, 0.1])
plt.ylim([0, 4])
plt.xlabel(feature)
plt.ylabel('-lnp')
plt.savefig('%s/volcano/%s_volcano_pw_rescale.pdf' % (output_dir, feature))
plt.close()

plt.subplots(figsize=(9, 6))
plt.scatter(df_mean[feature], df_p[feature], color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df_mean)):
    plt.text(x=df_mean[feature][i]+0.001, y=df_p[feature][i]+0.01, s=df_mean['gene'][i], size=5)
plt.scatter(df_mean_ctrl[feature], df_p_ctrl[feature], color=(0.85, 0.35, 0.25), s=20)
plt.xlim([-0.1, 0.1])
plt.ylim([0, 4])
plt.xlabel(feature)
plt.ylabel('-lnp')
plt.savefig('%s/volcano/%s_volcano_gene_rescale.pdf' % (output_dir, feature))
plt.close()

for feature in features2:
    plt.subplots(figsize=(9, 6))
    plt.scatter(df_mean_ex[feature], df_p_ex[feature], color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df_mean_ex)):
        plt.text(x=df_mean_ex[feature][i]+0.001, y=df_p_ex[feature][i]+0.01, s='%s(%s)' % (df_mean_ex['sample'][i], df_mean_ex['n'][i]), size=5)
    plt.scatter(df_mean_ctrl_ex[feature], df_p_ctrl_ex[feature], color=(0.85, 0.35, 0.25), s=20)
    # plt.axhline(y=-np.log(0.01), linestyle='--', color='red')
    limit = max(abs(df_mean_ex[feature].max()), abs(df_mean_ex[feature].min()))
    plt.xlim([-(limit+0.01), limit+0.01])
    plt.xlabel(feature)
    plt.ylabel('-lnp')
    plt.savefig('%s/volcano/%s_volcano_pw.pdf' % (output_dir, feature))
    plt.close()

for feature in features2:
    plt.subplots(figsize=(9, 6))
    plt.scatter(df_mean_ex[feature], df_p_ex[feature], color=(0.30, 0.30, 0.30), s=20)
    for i in range(len(df_mean_ex)):
        plt.text(x=df_mean_ex[feature][i]+0.001, y=df_p_ex[feature][i]+0.01, s='%s(%s)' % (df_mean_ex['gene'][i], df_mean_ex['n'][i]), size=5)
    plt.scatter(df_mean_ctrl_ex[feature], df_p_ctrl_ex[feature], color=(0.85, 0.35, 0.25), s=20)
    # plt.axhline(y=-np.log(0.01), linestyle='--', color='red')
    limit = max(abs(df_mean_ex[feature].max()), abs(df_mean_ex[feature].min()))
    plt.xlim([-(limit+0.01), limit+0.01])
    plt.xlabel(feature)
    plt.ylabel('-lnp')
    plt.savefig('%s/volcano/%s_volcano_gene.pdf' % (output_dir, feature))
    plt.close()