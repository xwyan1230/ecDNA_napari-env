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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230318_analysis_chemical-screen-nuclear_rep1_sp8_calibrate/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_nuclear_convex_dilation = -3
limit = 100
repeat = 50

# samples
exp = 'DM_6hr'
plates = ['DM_6hr']
rows = ['B', 'C', 'D', 'E', 'F', 'G']
columns = ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
num_total_samples = len(rows) * len(columns)

# load data
data = pd.DataFrame(columns=['exp', 'plate', 'well', 'p-w', 'compound', 'target', 'n_total', 'n_FOV',
                             'area_nuclear_mean', 'area_nuclear_mean_log2FC', 'area_nuclear_std', 'area_nuclear_p',
                             'mean_int_nuclear_mean', 'mean_int_nuclear_mean_log2FC', 'mean_int_nuclear_std', 'mean_int_nuclear_p',
                             'mean_int_MYC_mean', 'mean_int_MYC_mean_log2FC', 'mean_int_MYC_std', 'mean_int_MYC_p',
                             'mean_int_nuclear_cal_mean', 'mean_int_nuclear_cal_mean_log2FC', 'mean_int_nuclear_cal_std', 'mean_int_nuclear_cal_p',
                             'mean_int_MYC_cal_mean', 'mean_int_MYC_cal_mean_log2FC', 'mean_int_MYC_cal_std', 'mean_int_MYC_cal_p',
                             'total_int_nuclear_mean', 'total_int_nuclear_mean_log2FC', 'total_int_nuclear_std',
                             'total_int_nuclear_p',
                             'total_int_MYC_mean', 'total_int_MYC_mean_log2FC', 'total_int_MYC_std', 'total_int_MYC_p',
                             'total_int_nuclear_cal_mean', 'total_int_nuclear_cal_mean_log2FC',
                             'total_int_nuclear_cal_std', 'total_int_nuclear_cal_p',
                             'total_int_MYC_cal_mean', 'total_int_MYC_cal_mean_log2FC', 'total_int_MYC_cal_std',
                             'total_int_MYC_cal_p'])

df_WT = pd.read_csv('%s%s_n%s_WT.txt' % (data_dir1, exp, n_nuclear_convex_dilation), na_values=['.'], sep='\t')
df_WT['total_int_nuclear'] = df_WT['mean_int_nuclear'] * df_WT['area_nuclear']
df_WT['total_int_MYC'] = df_WT['mean_int_MYC'] * df_WT['area_nuclear']
df_WT['total_int_nuclear_cal'] = df_WT['mean_int_nuclear_cal'] * df_WT['area_nuclear']
df_WT['total_int_MYC_cal'] = df_WT['mean_int_MYC_cal'] * df_WT['area_nuclear']
df_gene = pd.read_csv('%sgene.txt' % data_dir1, na_values=['.'], sep='\t')


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
    df = pd.read_csv('%s%s_n%s.txt' % (data_dir1, plate, n_nuclear_convex_dilation), na_values=['.'], sep='\t')
    df['total_int_nuclear'] = df['mean_int_nuclear'] * df['area_nuclear']
    df['total_int_MYC'] = df['mean_int_MYC'] * df['area_nuclear']
    df['total_int_nuclear_cal'] = df['mean_int_nuclear_cal'] * df['area_nuclear']
    df['total_int_MYC_cal'] = df['mean_int_MYC_cal'] * df['area_nuclear']
    for j in range(num_total_samples):
        row = rows[int(j / len(columns))]
        column = columns[int(j - int(j / len(columns)) * len(columns))]
        print("%s: %s%s" % (plate, row, column))
        df_temp = df[df['well'] == '%s%s' % (row, column)].copy().reset_index(drop=True)
        area_nuclear_p = get_minuslnp(df_temp, df_WT, 'area_nuclear', limit, repeat)[1]
        mean_int_nuclear_p = get_minuslnp(df_temp, df_WT, 'mean_int_nuclear', limit, repeat)[1]
        mean_int_MYC_p = get_minuslnp(df_temp, df_WT, 'mean_int_MYC', limit, repeat)[1]
        mean_int_nuclear_cal_p = get_minuslnp(df_temp, df_WT, 'mean_int_nuclear_cal', limit, repeat)[1]
        mean_int_MYC_cal_p = get_minuslnp(df_temp, df_WT, 'mean_int_MYC_cal', limit, repeat)[1]
        total_int_nuclear_p = get_minuslnp(df_temp, df_WT, 'total_int_nuclear', limit, repeat)[1]
        total_int_MYC_p = get_minuslnp(df_temp, df_WT, 'total_int_MYC', limit, repeat)[1]
        total_int_nuclear_cal_p = get_minuslnp(df_temp, df_WT, 'total_int_nuclear_cal', limit, repeat)[1]
        total_int_MYC_cal_p = get_minuslnp(df_temp, df_WT, 'total_int_MYC_cal', limit, repeat)[1]
        df_gene_temp = df_gene[(df_gene['plate'] == i+1) & (df_gene['well'] == '%s%s' % (row, column))].copy().reset_index(drop=True)

        data.loc[len(data.index)] = [exp, plate, '%s%s' % (row, column), '%s-%s%s' % (i+1, row, column),
                                     df_gene_temp['compound'][0], df_gene_temp['target'][0],
                                     len(df_temp), len(df_temp)/max(df_temp['FOV']),
                                     np.mean(df_temp['area_nuclear']),
                                     np.log2(np.mean(df_temp['area_nuclear'])/np.mean(df_WT['area_nuclear'])),
                                     np.std(df_temp['area_nuclear']), area_nuclear_p,
                                     np.mean(df_temp['mean_int_nuclear']),
                                     np.log2(np.mean(df_temp['mean_int_nuclear']) / np.mean(df_WT['mean_int_nuclear'])),
                                     np.std(df_temp['mean_int_nuclear']),
                                     mean_int_nuclear_p,
                                     np.mean(df_temp['mean_int_MYC']),
                                     np.log2(np.mean(df_temp['mean_int_MYC']) / np.mean(df_WT['mean_int_MYC'])),
                                     np.std(df_temp['mean_int_MYC']), mean_int_MYC_p,
                                     np.mean(df_temp['mean_int_nuclear_cal']),
                                     np.log2(np.mean(df_temp['mean_int_nuclear_cal']) / np.mean(df_WT['mean_int_nuclear_cal'])),
                                     np.std(df_temp['mean_int_nuclear_cal']),
                                     mean_int_nuclear_cal_p,
                                     np.mean(df_temp['mean_int_MYC_cal']),
                                     np.log2(np.mean(df_temp['mean_int_MYC_cal']) / np.mean(df_WT['mean_int_MYC_cal'])),
                                     np.std(df_temp['mean_int_MYC_cal']),
                                     mean_int_MYC_cal_p,
                                     np.mean(df_temp['total_int_nuclear']),
                                     np.log2(np.mean(df_temp['total_int_nuclear']) / np.mean(df_WT['total_int_nuclear'])),
                                     np.std(df_temp['total_int_nuclear']),
                                     total_int_nuclear_p,
                                     np.mean(df_temp['total_int_MYC']),
                                     np.log2(np.mean(df_temp['total_int_MYC']) / np.mean(df_WT['total_int_MYC'])),
                                     np.std(df_temp['total_int_MYC']), total_int_MYC_p,
                                     np.mean(df_temp['total_int_nuclear_cal']),
                                     np.log2(np.mean(df_temp['total_int_nuclear_cal']) / np.mean(df_WT['total_int_nuclear_cal'])),
                                     np.std(df_temp['total_int_nuclear_cal']),
                                     total_int_nuclear_cal_p,
                                     np.mean(df_temp['total_int_MYC_cal']),
                                     np.log2(np.mean(df_temp['total_int_MYC_cal']) / np.mean(df_WT['total_int_MYC_cal'])),
                                     np.std(df_temp['total_int_MYC_cal']),
                                     total_int_MYC_cal_p]

data.to_csv('%s%s_n%s_sum.txt' % (output_dir, exp, n_nuclear_convex_dilation), index=False, sep='\t')
print("DONE!")

