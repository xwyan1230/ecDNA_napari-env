import napari
import shared.display as dis
import matplotlib.pyplot as plt
import seaborn as sns
import tifffile as tif
from skimage.measure import regionprops
import numpy as np
import seaborn as sns
import shared.segmentation as seg
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231104_analysis_ColoDMandHSR_Wee1_6hr/"

sample = 'ColoDM_Wee1_degrader_gH2AX_pH3'
conc_lst = ['XY06']
sample_name = 'DM_Wee1_degrader_-5'
export = 'N'

df = pd.DataFrame()
for conc in conc_lst:
    temp = pd.read_csv('%s/%s_%s.txt' % (master_folder, sample, conc), na_values=['.'], sep='\t')
    df = pd.concat([df, temp], axis=0)
df['ln_mean_int_gH2AX'] = np.log(df['gH2AX'])
df['ln_mean_int_green'] = np.log(df['green'])
df['ln_mean_int_red'] = np.log(df['red'])
df['ln_mean_int_pH3'] = np.log(df['pH3'])

df1 = df[df['ln_mean_int_pH3'] >= 6].copy().reset_index(drop=True)

red_cut = 6.1

print(len(df))
print(len(df1))
# mitotic
print("mitotic")
print(len(df1[df1['ln_mean_int_pH3']>8]))
print(len(df1[df1['ln_mean_int_pH3']>8])/len(df1))
# damaged
print("damaged")
print(len(df1[df1['ln_mean_int_gH2AX']>9]))
print(len(df1[df1['ln_mean_int_gH2AX']>9])/len(df1))
# double negative
print("negative")
print(len(df1[(df1['ln_mean_int_red'] < red_cut) & (df1['ln_mean_int_green'] < 6.5)]))
print(len(df1[(df1['ln_mean_int_red'] < red_cut) & (df1['ln_mean_int_green'] < 6.5)])/len(df1))
# S
print("S")
print(len(df1[(df1['ln_mean_int_red'] < red_cut) & (df1['ln_mean_int_green'] >= 6.5) & (df1['ln_mean_int_pH3'] <= 8)]))
print(len(df1[(df1['ln_mean_int_red'] < red_cut) & (df1['ln_mean_int_green'] >= 6.5) & (df1['ln_mean_int_pH3'] <= 8)])/len(df1))
# G1
print("G1")
print(len(df1[(df1['ln_mean_int_red'] >= red_cut) & (df1['ln_mean_int_green'] < 7) & (df1['ln_mean_int_pH3'] <= 8)]))
print(len(df1[(df1['ln_mean_int_red'] >= red_cut) & (df1['ln_mean_int_green'] < 7) & (df1['ln_mean_int_pH3'] <= 8)])/len(df1))
# G2
print("G2")
print(len(df1[(df1['ln_mean_int_red'] >= red_cut) & (df1['ln_mean_int_green'] >= 7) & (df1['ln_mean_int_pH3'] <= 8)]))
print(len(df1[(df1['ln_mean_int_red'] >= red_cut) & (df1['ln_mean_int_green'] >= 7) & (df1['ln_mean_int_pH3'] <= 8)])/len(df1))

"""hue_order = ['G1', 'S', 'G2', 'M']
line_colors = [(220/255, 20/255, 60/255), (154/255, 205/255, 50/255), (255/255, 127/255, 80/255), (199/255, 21/255, 133/255)]

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df1, y='ln_mean_int_green', x='ln_mean_int_red', s=20, alpha=0.5)
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df1, y='ln_mean_int_green', x='ln_mean_int_red', s=20)
sns.scatterplot(data=df1[(df1['ln_mean_int_red'] < red_cut) & (df1['ln_mean_int_green'] >= 6.5) & (df1['ln_mean_int_pH3'] <= 8)], y='ln_mean_int_green', x='ln_mean_int_red', s=20, color=(154/255, 205/255, 50/255))
sns.scatterplot(data=df1[(df1['ln_mean_int_red'] >= red_cut) & (df1['ln_mean_int_green'] < 7) & (df1['ln_mean_int_pH3'] <= 8)], y='ln_mean_int_green', x='ln_mean_int_red', s=20, color=(220/255, 20/255, 60/255))
sns.scatterplot(data=df1[(df1['ln_mean_int_red'] >= red_cut) & (df1['ln_mean_int_green'] >= 7) & (df1['ln_mean_int_pH3'] <= 8)], y='ln_mean_int_green', x='ln_mean_int_red', s=20, color=(255/255, 127/255, 80/255))
sns.scatterplot(data=df1[df1['ln_mean_int_pH3'] > 8], y='ln_mean_int_green', x='ln_mean_int_red', s=20, color=(199/255, 21/255, 133/255))
if export == 'Y':
    plt.savefig('%s/cellcycle_%s.pdf' % (master_folder, sample_name))
plt.show()


plt.subplots(figsize=(9, 6))
sns.distplot(df['ln_mean_int_pH3'], bins=50, kde=False, hist_kws={'range':(5,11)})
if export == 'Y':
    plt.savefig('%s/pH3_hist_%s.pdf' % (master_folder, sample_name))
plt.show()

plt.subplots(figsize=(9, 6))
sns.distplot(df['ln_mean_int_gH2AX'], bins=50, kde=False, hist_kws={'range':(7.5,11)})
if export == 'Y':
    plt.savefig('%s/gH2AX_hist_%s.pdf' % (master_folder, sample_name))
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df1, y='ln_mean_int_green', x='ln_mean_int_red', c=df1['ln_mean_int_pH3'].tolist(), cmap='coolwarm', s=20, alpha=1, vmin=5, vmax=11)
if export == 'Y':
    plt.savefig('%s/cellcycle_pH3_%s.pdf' % (master_folder, sample_name))
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df1, y='ln_mean_int_green', x='ln_mean_int_red', c=df1['ln_mean_int_gH2AX'].tolist(), cmap='coolwarm', s=20, alpha=1, vmin=8, vmax=10.5)
if export == 'Y':
    plt.savefig('%s/cellcycle_gH2AX_%s.pdf' % (master_folder, sample_name))
plt.show()"""


