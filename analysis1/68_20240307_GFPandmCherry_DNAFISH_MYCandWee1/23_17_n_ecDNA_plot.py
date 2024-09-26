import pandas as pd
import matplotlib.pyplot as plt
from shared.sinaplot import sinaplot
import numpy as np
import scipy.stats as stats
import seaborn as sns

# INPUT PARAMETERS
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"
sample = 'F8'

data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

df = pd.read_csv('%s/%s_group.txt' % (data_dir, sample), na_values=['.'], sep='\t')
df['ln_mean_int_DNAFISH'] = np.log(df['mean_int_DNAFISH'])

df_sort = df[df['group'].isin(['GFP', 'mCherry'])]

features = ['n_ecDNA', 'ln_mean_int_DNAFISH', 'total_area_ecDNA']

for feature in features:
    print(len(df_sort[df_sort['group'] == 'GFP']))
    print(len(df_sort[df_sort['group'] == 'mCherry']))
    print(np.mean(df_sort[df_sort['group'] == 'GFP'][feature]))
    print(np.mean(df_sort[df_sort['group'] == 'mCherry'][feature]))
    _, p_val = stats.ttest_ind(df_sort[df_sort['group'] == 'GFP'][feature], df_sort[df_sort['group'] == 'mCherry'][feature], equal_var=True)
    print(p_val)

    plt.subplots(figsize=(4, 9))
    sinaplot(data=df_sort, x='group', y=feature, alpha=0.7, violin=False, scale='area')
    plt.savefig('%s/%s/%s_%s_sinaplot.pdf' % (output_dir, sample, sample, feature))
    plt.show()

    plt.subplots(figsize=(4, 9))
    sns.violinplot(data=df_sort, x='group', y=feature, cut=0)
    plt.savefig('%s/%s/%s_%s_violin.pdf' % (output_dir, sample, sample, feature))
    plt.show()
