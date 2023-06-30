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
exp = 'HSR_6hr'
plates = ['HSR_6hr']
seq = pd.read_csv('%s%s_sequence.txt' % (data_dir, exp), na_values=['.'], sep='\t')
samples = seq['sample'].tolist()

# load data
data = pd.DataFrame(columns=['exp', 'plate', 'sample', 'well', 'p-w', 'compound', 'target', 'n_total', 'n_FOV',
                             'area_nuclear_mean', 'area_nuclear_std', 'area_nuclear_p',
                             'mean_int_nuclear_mean', 'mean_int_nuclear_std', 'mean_int_nuclear_p',
                             'mean_int_MYC_mean', 'mean_int_MYC_std', 'mean_int_MYC_p',
                             'mean_int_nuclear_cal_mean', 'mean_int_nuclear_cal_std', 'mean_int_nuclear_cal_p',
                             'mean_int_MYC_cal_mean', 'mean_int_MYC_cal_std', 'mean_int_MYC_cal_p'])

df_WT = pd.read_csv('%s%s_n%s_WT.txt' % (data_dir1, exp, n_nuclear_convex_dilation), na_values=['.'], sep='\t')
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
    for j in range(len(samples)):
        sample = samples[j]
        print("%s: %s" % (plate, sample))
        df_temp = df[df['sample'] == sample].copy().reset_index(drop=True)
        area_nuclear_p = get_minuslnp(df_temp, df_WT, 'area_nuclear', limit, repeat)[1]
        mean_int_nuclear_p = get_minuslnp(df_temp, df_WT, 'mean_int_nuclear', limit, repeat)[1]
        mean_int_MYC_p = get_minuslnp(df_temp, df_WT, 'mean_int_MYC', limit, repeat)[1]
        mean_int_nuclear_cal_p = get_minuslnp(df_temp, df_WT, 'mean_int_nuclear_cal', limit, repeat)[1]
        mean_int_MYC_cal_p = get_minuslnp(df_temp, df_WT, 'mean_int_MYC_cal', limit, repeat)[1]
        df_gene_temp = df_gene[(df_gene['plate'] == i+1) & (df_gene['well'] == sample.split('_')[0])].copy().reset_index(drop=True)

        data.loc[len(data.index)] = [exp, plate, sample, sample.split('_')[0], '%s-%s' % (i+1, sample),
                                     df_gene_temp['compound'][0], df_gene_temp['target'][0],
                                     len(df_temp), len(df_temp)/max(df_temp['FOV']),
                                     np.mean(df_temp['area_nuclear']), np.std(df_temp['area_nuclear']), area_nuclear_p,
                                     np.mean(df_temp['mean_int_nuclear']), np.std(df_temp['mean_int_nuclear']),
                                     mean_int_nuclear_p,
                                     np.mean(df_temp['mean_int_MYC']), np.std(df_temp['mean_int_MYC']), mean_int_MYC_p,
                                     np.mean(df_temp['mean_int_nuclear_cal']), np.std(df_temp['mean_int_nuclear_cal']),
                                     mean_int_nuclear_cal_p,
                                     np.mean(df_temp['mean_int_MYC_cal']), np.std(df_temp['mean_int_MYC_cal']),
                                     mean_int_MYC_cal_p]

data.to_csv('%s%s_n%s_sum_sample.txt' % (output_dir, exp, n_nuclear_convex_dilation), index=False, sep='\t')
print("DONE!")

