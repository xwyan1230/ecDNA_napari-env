import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.dataframe as dat
import shared.math as mat
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230712_analysis_chemical-screen_FUCCI/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

feature_lst = ['G1', 'S/G2', 'G2/M', 'G1/S']
name_lst = ['G1', 'SG2', 'G2M', 'G1S']
plate = 'HSR_FUCCI_2hr'

df = pd.read_csv('%s%s/summary.txt' % (data_dir, plate), na_values=['.'], sep='\t')
df_thresh = pd.read_csv('%sthresh.txt' % data_dir, na_values=['.'], sep='\t')
df_seq = pd.read_csv('%sseq.txt' % data_dir, na_values=['.'], sep='\t')

df['total'] = df['G1'] + df['S/G2'] + df['G2/M'] + df['G1/S']

df_ctrl = df[df['sample'].isin(['C3', 'C10', 'D6', 'F3', 'F10'])].copy().reset_index(drop=True)
G1_avg = sum(df_ctrl['G1'])/sum(df_ctrl['total'])
SG2_avg = sum(df_ctrl['S/G2'])/sum(df_ctrl['total'])
G2M_avg = sum(df_ctrl['G2/M'])/sum(df_ctrl['total'])
G1S_avg = sum(df_ctrl['G1/S'])/sum(df_ctrl['total'])
avg_lst = [G1_avg, SG2_avg, G2M_avg, G1S_avg]
print(G1_avg)
print(SG2_avg)
print(G2M_avg)
print(G1S_avg)

for f in range(len(feature_lst)):
    temp = []
    for i in range(len(df)):
        if df['total'][i] < 200:
            temp.append(0)
        else:
            phenotype = np.log2((df[feature_lst[f]][i]/df['total'][i]+0.00000001)/avg_lst[f])
            temp.append(phenotype)
    df['%s_log2FC' % feature_lst[f]] = temp
df.to_csv('%s%s/summary_log2FC.txt' % (output_dir, plate), index=False, sep='\t')

for f in range(len(feature_lst)):
    temp = []
    for i in range(len(df_seq)):
        temp.append(df[df['sample'] == df_seq['well'][i]]['%s_log2FC' % feature_lst[f]].tolist()[0])
    df_feature = pd.DataFrame(temp)
    df_feature.index = df_seq['compound'].tolist()
    df_feature.columns = [feature_lst[f]]

    fig, ax = plt.subplots(figsize=(9, 14))
    fig.subplots_adjust(left=0.3)
    ax1 = sns.heatmap(df_feature, cbar=0, linewidths=2, vmax=df_thresh[feature_lst[f]][1], vmin=df_thresh[feature_lst[f]][0], square=True, cmap='coolwarm', annot=False, fmt='.2f') # annot=True
    if not os.path.exists("%s%s/heatmap_line/" % (output_dir, plate)):
        os.makedirs("%s%s/heatmap_line/" % (output_dir, plate))
    plt.savefig('%s%s/heatmap_line/%s_seq.pdf' % (output_dir, plate, name_lst[f]))
    plt.close()
