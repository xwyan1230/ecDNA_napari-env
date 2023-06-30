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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230504_analysis_DM_2hr/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

# samples
ctrls = ['C3', 'C10', 'D6', 'F3', 'F10']
features = ['averageD', 'n_ecDNA']

# load data
df = pd.read_csv('%ssummary.txt' % data_dir1, na_values=['.'], sep='\t')
df_p = pd.read_csv('%ssummary_p.txt' % data_dir1, na_values=['.'], sep='\t')
df['averageD_p'] = df_p['averageD']
df['n_ecDNA_p'] = df_p['n_ecDNA']

df_WT = pd.DataFrame()
for i in range(len(ctrls)):
    df_temp = df[df['sample'] == ctrls[i]].copy().reset_index(drop=True)
    df_WT = pd.concat([df_WT, df_temp], axis=0).reset_index(drop=True)

# scatter plot

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='averageD_p', x='averageD_mean', color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df)):
    plt.text(x=df['averageD_mean'][i]+0.001, y=df['averageD_p'][i]+0.001, s=df['sample'][i], size=5)
sns.scatterplot(data=df_WT, y='averageD_p', x='averageD_mean', color=(0.85, 0.35, 0.25), s=20)
plt.savefig('%s/averageD_volcano_pw.pdf' % output_dir)
plt.close()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='n_ecDNA_p', x='n_ecDNA', color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df)):
    plt.text(x=df['n_ecDNA'][i]+0.1, y=df['n_ecDNA_p'][i]+0.1, s=df['sample'][i], size=5)
sns.scatterplot(data=df_WT, y='n_ecDNA_p', x='n_ecDNA', color=(0.85, 0.35, 0.25), s=20)
plt.savefig('%s/n_ecDNA_volcano_pw.pdf' % output_dir)
plt.close()

"""plt.subplots(figsize=(9, 6))
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
plt.close()"""