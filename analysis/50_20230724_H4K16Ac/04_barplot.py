import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.dataframe as dat
from shared.sinaplot import sinaplot
import shared.math as mat
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230724_analysis_H4K16Ac/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'HSR'
total_fov = 49
start_fov = 0
n = 0

hue_order = ['all', 'ecDNA', 'non_ecDNA']
line_colors = [(213/255, 213/255, 213/255), (199/255, 100/255, 106/255), (152/255, 179/255, 248/255),
               (195/255, 130/255, 130/255), (255/255, 100/255, 78/100)]
sns.set_palette(sns.color_palette(line_colors))

df = pd.read_csv(("%s%s_radial_new_n%s.txt" % (data_dir, sample, n)), na_values=['.'], sep='\t')

df_int = pd.DataFrame()
df_int['group'] = ['all'] * len(df) + ['ecDNA'] * len(df) + ['non_ecDNA'] * len(df)
df_int['int'] = df['mean_int_IF'].tolist() + df['mean_int_IF_ecDNA'].tolist() + df['mean_int_IF_nonecDNA'].tolist()

fig, ax = plt.subplots(figsize=(4, 9))
fig.subplots_adjust(left=0.2)
sns.barplot(data=df_int, x='group', y='int', order=hue_order)
plt.savefig('%s/%s_int.pdf' % (output_dir, sample))
plt.show()
