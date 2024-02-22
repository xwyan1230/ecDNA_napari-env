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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230514_analysis_mixing_Wee1-BRD4/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 4
limit = 250

# samples
exp = '20230512_mixing_Wee1-BRD4_6hr'
sample = '10_BRD4-GFP-6hr_Ctrl-mCh-6hr'
hue_order = ['GFP', 'mCherry']

data = pd.DataFrame(columns=['sample', 'limit',  'n',
                             'mean_area_nuclear', 'mean_total_area_ecDNA', 'mean_total_area_ratio_ecDNA', 'mean_mean_int_DNAFISH',
                             'mean_r1', 'mean_r2', 'mean_r3', 'mean_e1', 'mean_e2', 'mean_e3', 'mean_n_ecDNA',
                             'n_point2', 'mean_averageD_point2', 'std_averageD_point2', 'mean_n_ecDNA_point2', 'per_n_ecDNA10_point2'])

r1 = 23
r2 = 34
r3 = 37
e1 = 2
e2 = 5
e3 = 14


def get_minuslnp(data1: pd.DataFrame, data2: pd.DataFrame, feature: str, limit: int, repeat: int):
    minimum = np.min([len(data1), len(data2)])
    limit = minimum if (minimum < limit) else limit
    p_lst = []
    for j in range(repeat):
        p = -np.log(ks_2samp(data1[feature].sample(n=limit).tolist(), data2[feature].sample(n=limit).tolist())[1])
        p_lst.append(p)
    return limit, np.mean(p_lst)


df = pd.read_csv('%s/%s/%s_summary_new.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
feature = ['radial_curve_DNAFISH', 'radial_curve_edge_DNAFISH']
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
r1_GFP = [df_GFP['radial_curve_DNAFISH'][i][r1] for i in range(len(df_GFP))]
r1_mCherry = [df_mCherry['radial_curve_DNAFISH'][i][r1] for i in range(len(df_mCherry))]
r2_GFP = [df_GFP['radial_curve_DNAFISH'][i][r2] for i in range(len(df_GFP))]
r2_mCherry = [df_mCherry['radial_curve_DNAFISH'][i][r2] for i in range(len(df_mCherry))]
r3_GFP = [df_GFP['radial_curve_DNAFISH'][i][r3] for i in range(len(df_GFP))]
r3_mCherry = [df_mCherry['radial_curve_DNAFISH'][i][r3] for i in range(len(df_mCherry))]
e1_GFP = [df_GFP['radial_curve_edge_DNAFISH'][i][e1] for i in range(len(df_GFP))]
e1_mCherry = [df_mCherry['radial_curve_edge_DNAFISH'][i][e1] for i in range(len(df_mCherry))]
e2_GFP = [df_GFP['radial_curve_edge_DNAFISH'][i][e2] for i in range(len(df_GFP))]
e2_mCherry = [df_mCherry['radial_curve_edge_DNAFISH'][i][e2] for i in range(len(df_mCherry))]
e3_GFP = [df_GFP['radial_curve_edge_DNAFISH'][i][e3] for i in range(len(df_GFP))]
e3_mCherry = [df_mCherry['radial_curve_edge_DNAFISH'][i][e3] for i in range(len(df_mCherry))]
df_r_GFP = pd.DataFrame({'r1': r1_GFP, 'r2': r2_GFP, 'r3': r3_GFP, 'e1': e1_GFP, 'e2': e2_GFP, 'e3': e3_GFP})
df_r_mCherry = pd.DataFrame({'r1': r1_mCherry, 'r2': r2_mCherry, 'r3': r3_mCherry, 'e1': e1_mCherry, 'e2': e2_mCherry, 'e3': e3_mCherry})
mean_r1_GFP = np.mean(r1_GFP)
mean_r1_mCherry = np.mean(r1_mCherry)
mean_r2_GFP = np.mean(r2_GFP)
mean_r2_mCherry = np.mean(r2_mCherry)
mean_r3_GFP = np.mean(r3_GFP)
mean_r3_mCherry = np.mean(r3_mCherry)
mean_e1_GFP = np.mean(e1_GFP)
mean_e1_mCherry = np.mean(e1_mCherry)
mean_e2_GFP = np.mean(e2_GFP)
mean_e2_mCherry = np.mean(e2_mCherry)
mean_e3_GFP = np.mean(e3_GFP)
mean_e3_mCherry = np.mean(e3_mCherry)
p_r1 = get_minuslnp(df_r_mCherry, df_r_GFP, 'r1', limit, 50)[1]
p_r2 = get_minuslnp(df_r_mCherry, df_r_GFP, 'r2', limit, 50)[1]
p_r3 = get_minuslnp(df_r_mCherry, df_r_GFP, 'r3', limit, 50)[1]
p_e1 = get_minuslnp(df_r_mCherry, df_r_GFP, 'e1', limit, 50)[1]
p_e2 = get_minuslnp(df_r_mCherry, df_r_GFP, 'e2', limit, 50)[1]
p_e3 = get_minuslnp(df_r_mCherry, df_r_GFP, 'e3', limit, 50)[1]
mean_n_ecDNA_GFP = np.mean(df_GFP['n_ecDNA'])
mean_n_ecDNA_mCherry = np.mean(df_mCherry['n_ecDNA'])
p_n_ecDNA = get_minuslnp(df_GFP, df_mCherry, 'n_ecDNA', limit, 50)[1]

df_sort = df[df['total_area_ecDNA_sqrt_normalized'] > 0.2].copy().reset_index(drop=True)
df = df_sort
df_GFP = df[df['group'] == 'GFP'].copy().reset_index(drop=True)
df_mCherry = df[df['group'] == 'mCherry'].copy().reset_index(drop=True)

n_GFP_point2 = len(df_GFP)
n_mCherry_point2 = len(df_mCherry)

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
               mean_r1_mCherry, mean_r2_mCherry, mean_r3_mCherry, mean_e1_mCherry, mean_e2_mCherry, mean_e3_mCherry, mean_n_ecDNA_mCherry,
               n_mCherry_point2,
               mean_averageD_mCherry_point2, std_averageD_mCherry_point2,
               mean_n_ecDNA_mCherry_point2, per_n_ecDNA10_mCherry_point2]
data.loc[1] = ['%s_GFP' % sample, limit, n_GFP,
               mean_area_GFP, mean_total_area_ecDNA_GFP, mean_total_area_ratio_ecDNA_GFP, mean_mean_int_DNAFISH_GFP,
               mean_r1_GFP, mean_r2_GFP, mean_r3_GFP, mean_e1_GFP, mean_e2_GFP, mean_e3_GFP, mean_n_ecDNA_GFP,
               n_GFP_point2,
               mean_averageD_GFP_point2, std_averageD_GFP_point2,
               mean_n_ecDNA_GFP_point2, per_n_ecDNA10_GFP_point2]
data.loc[2] = ['p', limit, min([limit, n_GFP, n_mCherry]),
               p_area, p_total_area_ecDNA, p_total_area_ratio_ecDNA, p_mean_int_DNAFISH,
               p_r1, p_r2, p_r3, p_e1, p_e2, p_e3, p_n_ecDNA,
               min([limit, n_mCherry_point2, n_GFP_point2]),
               p_averageD_point2, p_averageD_point2,
               p_n_ecDNA_point2, p_n_ecDNA_point2]

data.to_csv('%s%s/%s_summary_forheatmap.txt' % (output_dir, sample, sample), index=False, sep='\t')
print("DONE!")
