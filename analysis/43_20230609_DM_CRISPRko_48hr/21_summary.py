import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
from shared.sinaplot import sinaplot
import shared.dataframe as dat
import seaborn as sns
from scipy.stats import ks_2samp
import numpy as np
import math
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 4
limit = 100

# samples
sample = 'E9'
hue_order = ['GFP', 'mCherry']

data = pd.DataFrame(columns=['sample', 'limit',  'n',
                             'mean_area_nuclear', 'mean_total_area_ecDNA', 'mean_total_area_ratio_ecDNA', 'mean_mean_int_DNAFISH',
                             'mean_r10', 'mean_r16', 'mean_n_ecDNA',
                             'n_point2', 'mean_r10_point2', 'mean_r16_point2',
                             'mean_averageD_point2', 'std_averageD_point2', 'mean_n_ecDNA_point2', 'per_n_ecDNA10_point2'])


def get_minuslnp(data1: pd.DataFrame, data2: pd.DataFrame, feature: str, limit: int, repeat: int):
    minimum = np.min([len(data1), len(data2)])
    limit = minimum if (minimum < limit) else limit
    p_lst = []
    for j in range(repeat):
        p = -np.log(ks_2samp(data1[feature].sample(n=limit).tolist(), data2[feature].sample(n=limit).tolist())[1])
        p_lst.append(p)
    return limit, np.mean(p_lst)


df = pd.read_csv('%s/%s/%s_n4_simplified.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
feature = ['radial_curve_DNAFISH']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]
df['r'] = np.sqrt(df['area_nuclear']/math.pi)
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA']/math.pi)
df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt']/df['r']
df['dis_to_hub_area_normalized'] = df['dis_to_hub_area']/df['r']
df_sample = df[df['group'].isin(hue_order)].copy().reset_index(drop=True)
df = df_sample

df_GFP = df[df['group'] == 'GFP'].copy().reset_index(drop=True)
df_mCherry = df[df['group'] == 'mCherry'].copy().reset_index(drop=True)

mean_area_GFP = np.mean(df_GFP['area_nuclear'])
mean_area_mCherry = np.mean(df_mCherry['area_nuclear'])
p_area = get_minuslnp(df_GFP, df_mCherry, 'area_nuclear', limit, 50)[1]
mean_total_area_ecDNA_GFP = np.mean(df_GFP['total_area_ecDNA'])
mean_total_area_ecDNA_mCherry = np.mean(df_mCherry['total_area_ecDNA'])
p_total_area_ecDNA = get_minuslnp(df_GFP, df_mCherry, 'total_area_ecDNA', limit, 50)[1]
mean_total_area_ratio_ecDNA_GFP = np.mean(df_GFP['total_area_ratio_ecDNA'])
mean_total_area_ratio_ecDNA_mCherry = np.mean(df_mCherry['total_area_ratio_ecDNA'])
p_total_area_ratio_ecDNA = get_minuslnp(df_GFP, df_mCherry, 'total_area_ratio_ecDNA', limit, 50)[1]
mean_mean_int_DNAFISH_GFP = np.mean(df_GFP['mean_int_DNAFISH'])
mean_mean_int_DNAFISH_mCherry = np.mean(df_mCherry['mean_int_DNAFISH'])
p_mean_int_DNAFISH = get_minuslnp(df_GFP, df_mCherry, 'mean_int_DNAFISH', limit, 50)[1]

n_GFP = len(df_GFP)
n_mCherry = len(df_mCherry)
r10_GFP = [df_GFP['radial_curve_DNAFISH'][i][10] for i in range(len(df_GFP))]
r10_mCherry = [df_mCherry['radial_curve_DNAFISH'][i][10] for i in range(len(df_mCherry))]
r16_GFP = [df_GFP['radial_curve_DNAFISH'][i][16] for i in range(len(df_GFP))]
r16_mCherry = [df_mCherry['radial_curve_DNAFISH'][i][16] for i in range(len(df_mCherry))]
df_r_GFP = pd.DataFrame({'r10': r10_GFP, 'r16': r16_GFP})
df_r_mCherry = pd.DataFrame({'r10': r10_mCherry, 'r16': r16_mCherry})
mean_r10_GFP = np.mean(r10_GFP)
mean_r10_mCherry = np.mean(r10_mCherry)
mean_r16_GFP = np.mean(r16_GFP)
mean_r16_mCherry = np.mean(r16_mCherry)
p_r10 = get_minuslnp(df_r_mCherry, df_r_GFP, 'r10', limit, 50)[1]
p_r16 = get_minuslnp(df_r_mCherry, df_r_GFP, 'r16', limit, 50)[1]
mean_n_ecDNA_GFP = np.mean(df_GFP['n_ecDNA'])
mean_n_ecDNA_mCherry = np.mean(df_mCherry['n_ecDNA'])
p_n_ecDNA = get_minuslnp(df_GFP, df_mCherry, 'n_ecDNA', limit, 50)[1]

df_sort = df[df['total_area_ecDNA_sqrt_normalized'] > 0.2].copy().reset_index(drop=True)
df = df_sort
df_GFP = df[df['group'] == 'GFP'].copy().reset_index(drop=True)
df_mCherry = df[df['group'] == 'mCherry'].copy().reset_index(drop=True)

n_GFP_point2 = len(df_GFP)
n_mCherry_point2 = len(df_mCherry)
r10_GFP_point2 = [df_GFP['radial_curve_DNAFISH'][i][10] for i in range(len(df_GFP))]
r10_mCherry_point2 = [df_mCherry['radial_curve_DNAFISH'][i][10] for i in range(len(df_mCherry))]
r16_GFP_point2 = [df_GFP['radial_curve_DNAFISH'][i][16] for i in range(len(df_GFP))]
r16_mCherry_point2 = [df_mCherry['radial_curve_DNAFISH'][i][16] for i in range(len(df_mCherry))]
df_r_GFP_point2 = pd.DataFrame({'r10': r10_GFP_point2, 'r16': r16_GFP_point2})
df_r_mCherry_point2 = pd.DataFrame({'r10': r10_mCherry_point2, 'r16': r16_mCherry_point2})
mean_r10_GFP_point2 = np.mean(r10_GFP_point2)
mean_r10_mCherry_point2 = np.mean(r10_mCherry_point2)
mean_r16_GFP_point2 = np.mean(r16_GFP_point2)
mean_r16_mCherry_point2 = np.mean(r16_mCherry_point2)
p_r10_point2 = get_minuslnp(df_r_mCherry_point2, df_r_GFP_point2, 'r10', limit, 50)[1]
p_r16_point2 = get_minuslnp(df_r_mCherry_point2, df_r_GFP_point2, 'r16', limit, 50)[1]

mean_averageD_GFP_point2 = np.mean(df_GFP['dis_to_hub_area_normalized'])
mean_averageD_mCherry_point2 = np.mean(df_mCherry['dis_to_hub_area_normalized'])
std_averageD_GFP_point2 = np.std(df_GFP['dis_to_hub_area_normalized'])
std_averageD_mCherry_point2 = np.std(df_mCherry['dis_to_hub_area_normalized'])
p_averageD_point2 = get_minuslnp(df_mCherry, df_GFP, 'dis_to_hub_area_normalized', limit, 50)[1]

mean_n_ecDNA_GFP_point2 = np.mean(df_GFP['n_ecDNA'])
mean_n_ecDNA_mCherry_point2 = np.mean(df_mCherry['n_ecDNA'])
per_n_ecDNA10_GFP_point2 = len(df_GFP[df_GFP['n_ecDNA'] > 10])/len(df_GFP)
per_n_ecDNA10_mCherry_point2 = len(df_mCherry[df_mCherry['n_ecDNA'] > 10])/len(df_mCherry)
p_n_ecDNA_point2 = get_minuslnp(df_mCherry, df_GFP, 'n_ecDNA', limit, 50)[1]

data.loc[0] = ['%s_mCherry' % sample, limit, n_mCherry,
               mean_area_mCherry, mean_total_area_ecDNA_mCherry, mean_total_area_ratio_ecDNA_mCherry, mean_mean_int_DNAFISH_mCherry,
               mean_r10_mCherry, mean_r16_mCherry, mean_n_ecDNA_mCherry,
               n_mCherry_point2, mean_r10_mCherry_point2, mean_r16_mCherry_point2,
               mean_averageD_mCherry_point2, std_averageD_mCherry_point2,
               mean_n_ecDNA_mCherry_point2, per_n_ecDNA10_mCherry_point2]
data.loc[1] = ['%s_GFP' % sample, limit, n_GFP,
               mean_area_GFP, mean_total_area_ecDNA_GFP, mean_total_area_ratio_ecDNA_GFP, mean_mean_int_DNAFISH_GFP,
               mean_r10_GFP, mean_r16_GFP, mean_n_ecDNA_GFP,
               n_GFP_point2, mean_r10_GFP_point2, mean_r16_GFP_point2,
               mean_averageD_GFP_point2, std_averageD_GFP_point2,
               mean_n_ecDNA_GFP_point2, per_n_ecDNA10_GFP_point2]
data.loc[2] = ['p', 50, min([limit, n_GFP, n_mCherry]),
               p_area, p_total_area_ecDNA, p_total_area_ratio_ecDNA, p_mean_int_DNAFISH,
               p_r10, p_r16, p_n_ecDNA,
               min([limit, n_mCherry_point2, n_GFP_point2]), p_r10_point2, p_r16_point2,
               p_averageD_point2, p_averageD_point2,
               p_n_ecDNA_point2, p_n_ecDNA_point2]

data.to_csv('%s%s/%s_summary.txt' % (output_dir, sample, sample), index=False, sep='\t')
print("DONE!")
