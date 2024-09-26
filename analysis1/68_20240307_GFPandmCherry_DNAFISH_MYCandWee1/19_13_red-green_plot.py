import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"
sample = 'F4'

data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

df = pd.read_csv('%s/%s.txt' % (data_dir, sample), na_values=['.'], sep='\t')
df['ln_mean_int_green'] = np.log(df['mean_int_green'])
df['ln_mean_int_red'] = np.log(df['mean_int_red'])

plt.subplots(figsize=(6, 6))
sns.scatterplot(data=df, y='ln_mean_int_green', x='ln_mean_int_red', color=(0.30, 0.30, 0.30), s=5, alpha=0.5)
if not os.path.exists("%s%s/" % (output_dir, sample)):
    os.makedirs("%s%s/" % (output_dir, sample))
plt.savefig('%s/%s/%s_red-green_scatter.pdf' % (output_dir, sample, sample))
plt.show()