import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
from shared.sinaplot import sinaplot
from scipy.stats import ks_2samp
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230313_analysis_chemical-screen-nuclear_rep1_EVOS/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

limit = 100
repeat = 50

# samples
exp = 'HSR_6hr'
plates = ['HSR_6hr']
rows = ['B', 'C', 'D', 'E', 'F', 'G']
columns = ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
num_total_samples = len(rows) * len(columns)

# load data
data = pd.DataFrame(columns=['exp', 'plate', 'well', 'p-w', 'compound', 'target', 'n_total', 'n_FOV',
                             'area_nuclear_mean', 'area_nuclear_std', 'area_nuclear_p',
                             'mean_int_nuclear_mean', 'mean_int_nuclear_std', 'mean_int_nuclear_p',
                             'mean_int_MYC_mean', 'mean_int_MYC_std', 'mean_int_MYC_p'])

df_WT = pd.read_csv('%s%s_WT.txt' % (data_dir, exp), na_values=['.'], sep='\t')
df_gene = pd.read_csv('%sgene.txt' % data_dir, na_values=['.'], sep='\t')


def get_minuslnp(data1: pd.DataFrame, data2: pd.DataFrame, feature: str, limit: int, repeat: int):
    minimum = np.min([len(data1), len(data2)])
    limit = minimum if (minimum < limit) else limit
    p_lst = []
    for j in range(repeat):
        p = -np.log(ks_2samp(data1[feature].sample(n=limit).tolist(), data2[feature].sample(n=limit).tolist())[1])
        p_lst.append(p)
    return limit, np.mean(p_lst)


for i in range(len(plates)):
    plate = plates[i]
    df = pd.read_csv('%s%s.txt' % (data_dir, plate), na_values=['.'], sep='\t')
    for j in range(num_total_samples):
        row = rows[int(j / len(columns))]
        column = columns[int(j - int(j / len(columns)) * len(columns))]
        print("%s: %s%s" % (plate, row, column))
        df_temp = df[df['well'] == '%s%s' % (row, column)].copy().reset_index(drop=True)
        area_nuclear_p = get_minuslnp(df_temp, df_WT, 'area_nuclear', limit, repeat)[1]
        mean_int_nuclear_p = get_minuslnp(df_temp, df_WT, 'mean_int_nuclear', limit, repeat)[1]
        mean_int_MYC_p = get_minuslnp(df_temp, df_WT, 'mean_int_MYC', limit, repeat)[1]
        df_gene_temp = df_gene[(df_gene['plate'] == i+1) & (df_gene['well'] == '%s%s' % (row, column))].copy().reset_index(drop=True)

        data.loc[len(data.index)] = [exp, plate, '%s%s' % (row, column), '%s-%s%s' % (i+1, row, column),
                                     df_gene_temp['compound'][0], df_gene_temp['target'][0],
                                     len(df_temp), len(df_temp)/max(df_temp['FOV']),
                                     np.mean(df_temp['area_nuclear']), np.std(df_temp['area_nuclear']), area_nuclear_p,
                                     np.mean(df_temp['mean_int_nuclear']), np.std(df_temp['mean_int_nuclear']),
                                     mean_int_nuclear_p,
                                     np.mean(df_temp['mean_int_MYC']), np.std(df_temp['mean_int_MYC']), mean_int_MYC_p]

data.to_csv('%s%s_sum.txt' % (output_dir, exp), index=False, sep='\t')
print("DONE!")

