import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 4

# samples
sample = 'C11'

df = pd.read_csv('%s/%s/%s_n4_simplified.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
df['ln_mean_int_green'] = np.log(df['mean_int_green'])
df['ln_mean_int_red'] = np.log(df['mean_int_red'])

df_align = pd.read_csv('%s/alignment/alignment.txt' % master_folder, na_values=['.'], sep='\t')
feature = ['GFP', 'mCherry']
for f in feature:
    df_align[f] = [dat.str_to_float(df_align[f][i]) for i in range(len(df_align))]

sample_lst = []
for i in range(len(df)):
    if (df['ln_mean_int_red'][i] < df_align[df_align['sample'] == sample]['GFP'].tolist()[0][0]) & (df['ln_mean_int_green'][i] > df_align[df_align['sample'] == sample]['GFP'].tolist()[0][1]):
        sample_lst.append('GFP')
    elif (df['ln_mean_int_red'][i] > df_align[df_align['sample'] == sample]['mCherry'].tolist()[0][0]) & (df['ln_mean_int_green'][i] < df_align[df_align['sample'] == sample]['mCherry'].tolist()[0][1]):
        sample_lst.append('mCherry')
    else:
        sample_lst.append('NA')
df['group'] = sample_lst

plt.subplots(figsize=(6, 6))
sns.scatterplot(data=df, y='ln_mean_int_green', x='ln_mean_int_red', color=(0.30, 0.30, 0.30), s=5, alpha=0.5)
sns.scatterplot(data=df[df['group'] == 'GFP'], y='ln_mean_int_green', x='ln_mean_int_red', color=(154/255, 205/255, 50/255), s=5, alpha=1)
sns.scatterplot(data=df[df['group'] == 'mCherry'], y='ln_mean_int_green', x='ln_mean_int_red', color=(220/255, 20/255, 60/255), s=5, alpha=1)
plt.savefig('%s/%s/%s_red-green_scatter_group.pdf' % (output_dir, sample, sample))
plt.show()

df.to_csv('%s%s/%s_n%s_simplified.txt' % (output_dir, sample, sample, n_dilation), index=False, sep='\t')
print("DONE!")