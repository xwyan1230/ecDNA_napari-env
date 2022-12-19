import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from skimage.measure import label, regionprops
import pandas as pd

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221215_analysis_radial-IF/20221117_immunoFISH_acid/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'H3K27me2me3'
n = 2
df_FISH = pd.read_csv(("%s%s/lineplot_csv/DNAFISH_lineplot%s.csv" % (data_dir, sample, n)), na_values=['.'], sep='\t')
df_IF = pd.read_csv(("%s%s/lineplot_csv/%s_lineplot%s.csv" % (data_dir, sample, sample, n)), na_values=['.'], sep='\t')

df_FISH['x'] = [df_FISH['X,Y'][i].split(',')[0] for i in range(len(df_FISH))]
df_FISH['y'] = [df_FISH['X,Y'][i].split(',')[1] for i in range(len(df_FISH))]
df_FISH['x'] = df_FISH['x'].astype(float)
df_FISH['y'] = df_FISH['y'].astype(float)
df_IF['x'] = [df_IF['X,Y'][i].split(',')[0] for i in range(len(df_IF))]
df_IF['y'] = [df_IF['X,Y'][i].split(',')[1] for i in range(len(df_IF))]
df_IF['x'] = df_IF['x'].astype(float)
df_IF['y'] = df_IF['y'].astype(float)

plt.subplots(figsize=(9, 6))
plt.plot(df_FISH['x'], df_FISH['y'], c='green', label='ecDNA')
plt.plot(df_IF['x'], df_IF['y'], c='red', label=sample)
plt.xlabel('Distance')
plt.ylabel('Intensity')
plt.legend()
plt.savefig('%slineplot_%s_%s.pdf' % (output_dir, sample, n))
plt.show()