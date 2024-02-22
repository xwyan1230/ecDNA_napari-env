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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221215_analysis_radial-IF/20221117_immunoFISH_acid/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'H3K27me2me3'
n = 0
hue_order = [2, 1]

df = pd.read_csv(("%s%s_radial_new_n%s.txt" % (data_dir, sample, n)), na_values=['.'], sep='\t')

feature_lst = ['IF_int', 'DNAFISH_seg_label']
for f in feature_lst:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

df_int = pd.DataFrame()
seg_lst = []
int_lst = []
for i in range(len(df)):
    seg_lst = seg_lst + df['DNAFISH_seg_label'][i] + [2] * len(df['DNAFISH_seg_label'][i])
    int_lst = int_lst + df['IF_int'][i] + df['IF_int'][i]
df_int['seg'] = seg_lst
df_int['int'] = int_lst
df_int.to_csv('%s%s_int.txt' % (output_dir, sample), index=False, sep='\t')
print("DONE!")

"""fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
sinaplot(data=df_int, x='seg', y='int', order=hue_order, violin=False, scale='area', point_size=2)
plt.savefig('%s/%s_int.pdf' % (output_dir, sample))
plt.show()"""
