import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
from shared.sinaplot import sinaplot
import seaborn as sns
import numpy as np
import math
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230514_analysis_mixing_Wee1-BRD4/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 4

# samples
df_seq = pd.read_csv('%sseq_mCh.txt' % data_dir, na_values=['.'], sep='\t')

col = ['sample', 'n', 'mean_area_nuclear', 'mean_total_area_ecDNA', 'mean_total_area_ratio_ecDNA',
       'mean_mean_int_DNAFISH', 'mean_r1', 'mean_r2', 'mean_r3', 'mean_e1', 'mean_e2', 'mean_e3', 'mean_n_ecDNA',
       'n_point2',
       'mean_averageD_point2', 'std_averageD_point2', 'mean_n_ecDNA_point2', 'per_n_ecDNA10_point2', 'gene']
data_mean = pd.DataFrame(columns=col)
data_p = pd.DataFrame(columns=col)


def get_phenotype(df, feature):
    return np.log2(df[feature][0]/(df[feature][1]+0.0001))

i = 0
j = 0
for k in range(len(df_seq)):
    sample = df_seq['location'][k]
    df = pd.read_csv('%s/%s/%s_summary_forheatmap.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
    data_mean.loc[len(data_mean.index)] = [sample, df['n'][2], get_phenotype(df, 'mean_area_nuclear'),
                                           get_phenotype(df, 'mean_total_area_ecDNA'),
                                           get_phenotype(df, 'mean_total_area_ratio_ecDNA'),
                                           get_phenotype(df, 'mean_mean_int_DNAFISH'),
                                           get_phenotype(df, 'mean_r1'),
                                           get_phenotype(df, 'mean_r2'),
                                           get_phenotype(df, 'mean_r3'),
                                           get_phenotype(df, 'mean_e1'),
                                           get_phenotype(df, 'mean_e2'),
                                           get_phenotype(df, 'mean_e3'),
                                           get_phenotype(df, 'mean_n_ecDNA'),
                                           df['n_point2'][2],
                                           get_phenotype(df, 'mean_averageD_point2'),
                                           get_phenotype(df, 'std_averageD_point2'),
                                           get_phenotype(df, 'mean_n_ecDNA_point2'),
                                           get_phenotype(df, 'per_n_ecDNA10_point2'),
                                           df_seq['gene'][k]]
    data_p.loc[len(data_p.index)] = [sample, df['n'][2], df['mean_area_nuclear'][2],
                                           df['mean_total_area_ecDNA'][2],
                                           df['mean_total_area_ratio_ecDNA'][2],
                                           df['mean_mean_int_DNAFISH'][2],
                                           df['mean_r1'][2],
                                           df['mean_r2'][2],
                                           df['mean_r3'][2],
                                           df['mean_e1'][2],
                                     df['mean_e2'][2], df['mean_e3'][2],
                                           df['mean_n_ecDNA'][2],
                                           df['n_point2'][2],
                                           df['mean_averageD_point2'][2],
                                           df['std_averageD_point2'][2],
                                           df['mean_n_ecDNA_point2'][2],
                                           df['per_n_ecDNA10_point2'][2],
                                           df_seq['gene'][k]]

data_mean.to_csv('%ssummary_mean_mCh.txt' % output_dir, index=False, sep='\t')
data_p.to_csv('%ssummary_p_mCh.txt' % output_dir, index=False, sep='\t')

# GFP
df_seq = pd.read_csv('%sseq_GFP.txt' % data_dir, na_values=['.'], sep='\t')

col = ['sample', 'n', 'mean_area_nuclear', 'mean_total_area_ecDNA', 'mean_total_area_ratio_ecDNA',
       'mean_mean_int_DNAFISH', 'mean_r1', 'mean_r2', 'mean_r3', 'mean_e1', 'mean_e2', 'mean_e3', 'mean_n_ecDNA',
       'n_point2',
       'mean_averageD_point2', 'std_averageD_point2', 'mean_n_ecDNA_point2', 'per_n_ecDNA10_point2', 'gene']
data_mean = pd.DataFrame(columns=col)
data_p = pd.DataFrame(columns=col)


def get_phenotype(df, feature):
    return np.log2(df[feature][1]/(df[feature][0]+0.0001))

i = 0
j = 0
for k in range(len(df_seq)):
    sample = df_seq['location'][k]
    df = pd.read_csv('%s/%s/%s_summary_forheatmap.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
    data_mean.loc[len(data_mean.index)] = [sample, df['n'][2], get_phenotype(df, 'mean_area_nuclear'),
                                           get_phenotype(df, 'mean_total_area_ecDNA'),
                                           get_phenotype(df, 'mean_total_area_ratio_ecDNA'),
                                           get_phenotype(df, 'mean_mean_int_DNAFISH'),
                                           get_phenotype(df, 'mean_r1'),
                                           get_phenotype(df, 'mean_r2'),
                                           get_phenotype(df, 'mean_r3'),
                                           get_phenotype(df, 'mean_e1'),
                                           get_phenotype(df, 'mean_e2'),
                                           get_phenotype(df, 'mean_e3'),
                                           get_phenotype(df, 'mean_n_ecDNA'),
                                           df['n_point2'][2],
                                           get_phenotype(df, 'mean_averageD_point2'),
                                           get_phenotype(df, 'std_averageD_point2'),
                                           get_phenotype(df, 'mean_n_ecDNA_point2'),
                                           get_phenotype(df, 'per_n_ecDNA10_point2'),
                                           df_seq['gene'][k]]
    data_p.loc[len(data_p.index)] = [sample, df['n'][2], df['mean_area_nuclear'][2],
                                           df['mean_total_area_ecDNA'][2],
                                           df['mean_total_area_ratio_ecDNA'][2],
                                           df['mean_mean_int_DNAFISH'][2],
                                           df['mean_r1'][2],
                                           df['mean_r2'][2],
                                           df['mean_r3'][2],
                                           df['mean_e1'][2],
                                     df['mean_e2'][2], df['mean_e3'][2],
                                           df['mean_n_ecDNA'][2],
                                           df['n_point2'][2],
                                           df['mean_averageD_point2'][2],
                                           df['std_averageD_point2'][2],
                                           df['mean_n_ecDNA_point2'][2],
                                           df['per_n_ecDNA10_point2'][2],
                                           df_seq['gene'][k]]

data_mean.to_csv('%ssummary_mean_GFP.txt' % output_dir, index=False, sep='\t')
data_p.to_csv('%ssummary_p_GFP.txt' % output_dir, index=False, sep='\t')

print("DONE!")
