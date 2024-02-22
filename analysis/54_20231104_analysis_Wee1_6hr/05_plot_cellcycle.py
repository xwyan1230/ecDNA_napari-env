import napari
import shared.display as dis
import matplotlib.pyplot as plt
import seaborn as sns
import tifffile as tif
from skimage.measure import regionprops
import numpy as np
import seaborn as sns
from shared.sinaplot import sinaplot
import shared.segmentation as seg
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231104_analysis_ColoDMandHSR_Wee1_6hr/"

sample = 'ColoDM_Wee1_degrader_gH2AX_pH3'
conc_lst = [['XY05', 'XY06'],
            ['XY10', 'XY11', 'XY12', 'XY13', 'XY14'],
            ['XY15', 'XY16', 'XY17'],
            ['XY23', 'XY24', 'XY25'],
            ['XY29', 'XY30', 'XY31', 'XY33', 'XY35', 'XY37'],
            ['XY38', 'XY39', 'XY43', 'XY44'],
            ['XY46', 'XY47', 'XY49', 'XY50', 'XY51', 'XY52'],
            ['XY53', 'XY54', 'XY55', 'XY58', 'XY59', 'XY60'],
            ['XY61', 'XY62', 'XY63'],
            ['XY68']]
"""conc_lst = [['XY05'],
            ['XY09', 'XY10'],
            ['XY13'],
            ['XY17', 'XY18', 'XY19', 'XY20'],
            ['XY21', 'XY22'],
            ['XY25', 'XY26', 'XY27', 'XY28'],
            ['XY29', 'XY30', 'XY31', 'XY32'],
            ['XY33', 'XY34', 'XY35', 'XY36'],
            ['XY37', 'XY38'],
            ['XY41', 'XY42']]"""
group_name = ['-5', '-4', '-3', '-2', '-1', '0', '1', '2', '3', '4']
sample_name = 'DM_Wee1_degrader'

line_colors = [(154/255, 205/255, 50/255)] * 10
# line_colors = [(255/255, 127/255, 80/255)] * 10
sns.set_palette(sns.color_palette(line_colors))

df = pd.DataFrame()
for n_group in range(len(conc_lst)):
    group_df = pd.DataFrame()
    for conc in conc_lst[n_group]:
        temp = pd.read_csv('%s/%s_%s_cellcycle.txt' % (master_folder, sample, conc), na_values=['.'], sep='\t')
        group_df = pd.concat([group_df, temp], axis=0)
    group_df['group'] = [group_name[n_group]]*len(group_df)
    df = pd.concat([df, group_df], axis=0)

df2 = df[df['ln_mean_int_pH3'] >= 6].copy().reset_index(drop=True)
df1 = df[(df['ln_mean_int_pH3'] >= 6) & (df['cellcycle'] == 'M')].copy().reset_index(drop=True)

pos_lst = []
for i in group_name:
    pos_lst.append(len(df1[(df1['group'] == i) & (df1['ln_mean_int_gH2AX'] > 9)])/ len(df2[df2['group'] == i]))
print(pos_lst)

"""plt.subplots(figsize=(9, 6))
sns.barplot(y=pos_lst, x=group_name)
plt.ylim([0, 1.0])
plt.savefig('%s/gH2AX_bar_M_%s.pdf' % (master_folder, sample_name))
plt.show()"""


"""fig, ax = plt.subplots(figsize=(15, 6))
fig.subplots_adjust(left=0.2)
feature = 'ln_mean_int_gH2AX'
sinaplot(data=df1, x='group', y=feature, order=group_name, violin=False, scale='area', point_size=2)
plt.ylim([7.5, 10.5])
plt.savefig('%s/gH2AX_M_%s.pdf' % (master_folder, sample_name))
plt.show()"""