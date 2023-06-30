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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 4

# samples
rows = ['C', 'D']
columns = ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
additional_sample = ['E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9']
total = len(rows) * len(columns) + len(additional_sample)

col = ['sample', 'n', 'mean_area_nuclear', 'mean_total_area_ecDNA', 'mean_total_area_ratio_ecDNA',
       'mean_mean_int_DNAFISH', 'mean_r10', 'mean_r16', 'mean_n_ecDNA', 'n_point2', 'mean_r10_point2', 'mean_r16_point2',
       'mean_averageD_point2', 'std_averageD_point2', 'mean_n_ecDNA_point2', 'per_n_ecDNA10_point2', 'gene']
data_mean = pd.DataFrame(columns=col)
data_p = pd.DataFrame(columns=col)

gene = pd.read_csv('%sgenelist.txt' % data_dir, na_values=['.'], sep='\t')


def get_phenotype(df, feature):
    return np.log2(df[feature][0]/(df[feature][1]+0.0001))

i = 0
j = 0
for k in range(total):
    if len(columns)*i + j+1 < (len(rows) * len(columns)):
        i = int(k / len(columns))
        j = int(k % len(columns))
        sample = '%s%s' % (rows[i], columns[j])
    else:
        sample = additional_sample[k-len(rows) * len(columns)]
    df = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
    data_mean.loc[len(data_mean.index)] = [sample, df['n'][2], get_phenotype(df, 'mean_area_nuclear'),
                                           get_phenotype(df, 'mean_total_area_ecDNA'),
                                           get_phenotype(df, 'mean_total_area_ratio_ecDNA'),
                                           get_phenotype(df, 'mean_mean_int_DNAFISH'),
                                           get_phenotype(df, 'mean_r10'),
                                           get_phenotype(df, 'mean_r16'),
                                           get_phenotype(df, 'mean_n_ecDNA'),
                                           df['n_point2'][2], get_phenotype(df, 'mean_r10_point2'),
                                           get_phenotype(df, 'mean_r16_point2'),
                                           get_phenotype(df, 'mean_averageD_point2'),
                                           get_phenotype(df, 'std_averageD_point2'),
                                           get_phenotype(df, 'mean_n_ecDNA_point2'),
                                           get_phenotype(df, 'per_n_ecDNA10_point2'),
                                           gene[gene['location'] == sample]['gene'].tolist()[0]]
    data_p.loc[len(data_p.index)] = [sample, df['n'][2], df['mean_area_nuclear'][2],
                                           df['mean_total_area_ecDNA'][2],
                                           df['mean_total_area_ratio_ecDNA'][2],
                                           df['mean_mean_int_DNAFISH'][2],
                                           df['mean_r10'][2],
                                           df['mean_r16'][2],
                                           df['mean_n_ecDNA'][2],
                                           df['n_point2'][2], df['mean_r10_point2'][2],
                                           df['mean_r16_point2'][2],
                                           df['mean_averageD_point2'][2],
                                           df['std_averageD_point2'][2],
                                           df['mean_n_ecDNA_point2'][2],
                                           df['per_n_ecDNA10_point2'][2],
                                           gene[gene['location'] == sample]['gene'].tolist()[0]]

data_mean.to_csv('%ssummary_mean.txt' % output_dir, index=False, sep='\t')
data_p.to_csv('%ssummary_p.txt' % output_dir, index=False, sep='\t')

print("DONE!")
