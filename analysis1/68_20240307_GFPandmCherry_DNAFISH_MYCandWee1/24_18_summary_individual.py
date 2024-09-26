import pandas as pd
from scipy.stats import ks_2samp
import numpy as np

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"
sample = 'F8'

data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

limit = 100

# samples
hue_order = ['GFP', 'mCherry']

data = pd.DataFrame(columns=['sample', 'limit',  'n',
                             'mean_area_nuclear', 'mean_total_area_ecDNA', 'mean_total_area_ratio_ecDNA',
                             'mean_mean_int_DNAFISH', 'mean_n_ecDNA', 'peak_relativer', 'fwhm_relativer',
                             'center_relativer', 'peak_value_relativer', 'peak_absoluter', 'fwhm_absoluter',
                             'center_absoluter', 'peak_value_absoluter'])


def get_minuslnp(data1: pd.DataFrame, data2: pd.DataFrame, feature: str, limit: int, repeat: int):
    minimum = np.min([len(data1), len(data2)])
    limit = minimum if (minimum < limit) else limit
    p_lst = []
    for j in range(repeat):
        p = -np.log(ks_2samp(data1[feature].sample(n=limit).tolist(), data2[feature].sample(n=limit).tolist())[1])
        p_lst.append(p)
    return limit, np.mean(p_lst)


def get_max_index(lst: list):
    return lst.index(max(lst))


def get_surrounding_points(lst: list, y):
    out = []
    for i in range(len(lst)):
        if (y < lst[i]) & (len(out) == 0):
            out.append([i-1, lst[i-1]])
            out.append([i, lst[i]])
        elif (y > lst[i]) & (len(out) == 2):
            out.append([i - 1, lst[i - 1]])
            out.append([i, lst[i]])
    return out


def find_x_from_y_and_line(point1: list, point2: list, y):
    return (y-point1[1])*(point2[0]-point1[0])/(point2[1]-point1[1]) + point1[0]


def get_fwhm_and_center(lst: list):
    halfpeak_value = 0.5*max(lst)
    surrounding_points = get_surrounding_points(lst, halfpeak_value)
    x1 = find_x_from_y_and_line(surrounding_points[0], surrounding_points[1], halfpeak_value)
    x2 = find_x_from_y_and_line(surrounding_points[2], surrounding_points[3], halfpeak_value)
    print(x1)
    print(x2)
    fwhm = x2-x1
    center = (x2+x1)/2
    print(fwhm)
    print(center)
    return fwhm, center


df = pd.read_csv('%s/%s_group.txt' % (data_dir, sample), na_values=['.'], sep='\t')
df_relativer = pd.read_csv('%s/%s/%s_radial_summary_relativer.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
df_absoluter = pd.read_csv('%s/%s/%s_radial_summary_absoluter.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')

mCherry_relativer_lst = df_relativer.iloc[1].tolist()[20:-1] + [0]
GFP_relativer_lst = df_relativer.iloc[0].tolist()[20:-1] + [0]
mCherry_absoluter_lst = df_absoluter.iloc[1].tolist()[20:-1] + [0]
GFP_absoluter_lst = df_absoluter.iloc[0].tolist()[20:-1] + [0]

mCherry_relativer_peak_value = max(mCherry_relativer_lst)
GFP_relativer_peak_value = max(GFP_relativer_lst)
mCherry_absoluter_peak_value = max(mCherry_absoluter_lst)
GFP_absoluter_peak_value = max(GFP_absoluter_lst)

mCherry_relativer_peak = get_max_index(mCherry_relativer_lst)+20
GFP_relativer_peak = get_max_index(GFP_relativer_lst)+20
mCherry_absoluter_peak = get_max_index(mCherry_absoluter_lst)+20
GFP_absoluter_peak = get_max_index(GFP_absoluter_lst)+20

mCherry_relativer_fwhm, mCherry_relativer_center = get_fwhm_and_center(mCherry_relativer_lst)
GFP_relativer_fwhm, GFP_relativer_center = get_fwhm_and_center(GFP_relativer_lst)
mCherry_absoluter_fwhm, mCherry_absoluter_center = get_fwhm_and_center(mCherry_absoluter_lst)
GFP_absoluter_fwhm, GFP_absoluter_center = get_fwhm_and_center(GFP_absoluter_lst)

mCherry_relativer_center = mCherry_relativer_center + 20
GFP_relativer_center = GFP_relativer_center + 20
mCherry_absoluter_center = mCherry_absoluter_center + 20
GFP_absoluter_center = GFP_absoluter_center + 20

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
mean_n_ecDNA_GFP = np.mean(df_GFP['n_ecDNA'])
mean_n_ecDNA_mCherry = np.mean(df_mCherry['n_ecDNA'])
p_n_ecDNA = get_minuslnp(df_GFP, df_mCherry, 'n_ecDNA', limit, 50)[1]

data.loc[0] = ['%s_mCherry' % sample, limit, n_mCherry,
               mean_area_mCherry, mean_total_area_ecDNA_mCherry, mean_total_area_ratio_ecDNA_mCherry, mean_mean_int_DNAFISH_mCherry,
               mean_n_ecDNA_mCherry, mCherry_relativer_peak, mCherry_relativer_fwhm, mCherry_relativer_center, mCherry_relativer_peak_value,
               mCherry_absoluter_peak, mCherry_absoluter_fwhm, mCherry_absoluter_center, mCherry_absoluter_peak_value]
data.loc[1] = ['%s_GFP' % sample, limit, n_GFP,
               mean_area_GFP, mean_total_area_ecDNA_GFP, mean_total_area_ratio_ecDNA_GFP, mean_mean_int_DNAFISH_GFP,
               mean_n_ecDNA_GFP, GFP_relativer_peak, GFP_relativer_fwhm, GFP_relativer_center, GFP_relativer_peak_value,
               GFP_absoluter_peak, GFP_absoluter_fwhm, GFP_absoluter_center, GFP_absoluter_peak_value]
data.loc[2] = ['p', limit, min([limit, n_GFP, n_mCherry]),
               p_area, p_total_area_ecDNA, p_total_area_ratio_ecDNA, p_mean_int_DNAFISH,
               p_n_ecDNA, 0, 0, 0, 0, 0, 0, 0, 0]

data.to_csv('%s%s_summary.txt' % (output_dir, sample), index=False, sep='\t')
print("DONE!")
