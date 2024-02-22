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
group = 'mCherry'
row = 1

# samples
rows = ['C', 'D', 'E', 'F', 'G']
columns = ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
additional_sample = []
total = len(rows) * len(columns) + len(additional_sample)

col = ['sample', 'n', 'mean_area_nuclear', 'mean_total_area_ecDNA', 'mean_total_area_ratio_ecDNA',
       'mean_mean_int_DNAFISH', 'peak_relativer', 'fwhm_relativer', 'center_relativer', 'peak_absoluter',
       'fwhm_absoluter', 'center_absoluter', 'mean_n_ecDNA',
       'n_point2',
       'mean_averageD_point2', 'std_averageD_point2', 'mean_n_ecDNA_point2', 'per_n_ecDNA10_point2', 'gene']

data = pd.DataFrame(columns=col)

gene = pd.read_csv('%sgenelist.txt' % data_dir, na_values=['.'], sep='\t')


def get_phenotype_log2(df, feature):
    return np.log2(df[feature][0]/(df[feature][1]+0.0001))


def get_phenotype(df, feature):
    return df[feature][0]-df[feature][1]


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
    data.loc[len(data.index)] = [sample, df['n'][row], df['mean_area_nuclear'][row],
                                 df['mean_total_area_ecDNA'][row],
                                 df['mean_total_area_ratio_ecDNA'][row],
                                 df['mean_mean_int_DNAFISH'][row],
                                 df['peak_relativer'][row],
                                 df['fwhm_relativer'][row],
                                 df['center_relativer'][row],
                                 df['peak_absoluter'][row],
                                 df['fwhm_absoluter'][row], df['center_absoluter'][row],
                                 df['mean_n_ecDNA'][row],
                                 df['n_point2'][row],
                                 df['mean_averageD_point2'][row],
                                 df['std_averageD_point2'][row],
                                 df['mean_n_ecDNA_point2'][row],
                                 df['per_n_ecDNA10_point2'][row],
                                 gene[gene['location'] == sample]['gene'].tolist()[0]]

data.to_csv('%ssummary_%s.txt' % (output_dir, group), index=False, sep='\t')

print("DONE!")
