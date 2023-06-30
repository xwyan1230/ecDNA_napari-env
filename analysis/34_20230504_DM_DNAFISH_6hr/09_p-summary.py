import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.dataframe as dat
from scipy.stats import ks_2samp
import shared.math as mat
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230504_analysis_DM_6hr/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

feature_lst = ['sample', 'averageD', 'n_ecDNA']
row_lst = ['B', 'C', 'D', 'E', 'F', 'G']
column_lst = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
sample_lst = []
for row in row_lst:
    sample_lst = sample_lst + ['%s%s' % (row, column) for column in column_lst]
ctrl_lst = ['C3', 'C10', 'D6', 'F3', 'F10']
df_ctrl = pd.DataFrame()
for ctrl in ctrl_lst:
    df = pd.read_csv('%s/txt/%s_n4.txt' % (data_dir1, ctrl), na_values=['.'], sep='\t')
    if len(df) >= 100:
        df_ctrl = pd.concat([df_ctrl, df.sample(n=100).reset_index()], axis=0)
    else:
        df_ctrl = pd.concat([df_ctrl, df.sample(n=len(df)).reset_index()], axis=0)
df_ctrl['r'] = np.sqrt(df_ctrl['area_nuclear'] / math.pi)
df_ctrl['total_area_ecDNA_sqrt'] = np.sqrt(df_ctrl['total_area_ecDNA'] / math.pi)
df_ctrl['total_area_ecDNA_sqrt_normalized'] = df_ctrl['total_area_ecDNA_sqrt'] / df_ctrl['r']
df_ctrl['dis_to_hub_area_normalized'] = df_ctrl['dis_to_hub_area'] / df_ctrl['r']


def get_minuslnp(data1: pd.DataFrame, data2: pd.DataFrame, feature: str, limit: int, repeat: int):
    minimum = np.min([len(data1), len(data2)])
    limit = minimum if (minimum < limit) else limit
    p_lst = []
    for j in range(repeat):
        p = -np.log(ks_2samp(data1[feature].sample(n=limit).tolist(), data2[feature].sample(n=limit).tolist())[1])
        p_lst.append(p)
    return limit, np.mean(p_lst)


df_p = pd.DataFrame(columns=feature_lst)
for sample in sample_lst:
    if os.path.exists('%s/txt/%s_n4.txt' % (data_dir1, sample)):
        df = pd.read_csv('%s/txt/%s_n4.txt' % (data_dir1, sample), na_values=['.'], sep='\t')
        df['r'] = np.sqrt(df['area_nuclear'] / math.pi)
        df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA'] / math.pi)
        df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt'] / df['r']
        df['dis_to_hub_area_normalized'] = df['dis_to_hub_area'] / df['r']
        averageD_p = get_minuslnp(df, df_ctrl, 'dis_to_hub_area_normalized', 100, 50)[1]
        n_ecDNA_p = get_minuslnp(df, df_ctrl, 'n_ecDNA', 100, 50)[1]
        df_p.loc[len(df_p.index)] = [sample, averageD_p, n_ecDNA_p]
    else:
        df_p.loc[len(df_p.index)] = [sample, 0, 0]
df_p.to_csv('%s/summary_p.txt' % output_dir, index=False, sep='\t')
print("DONE!")
