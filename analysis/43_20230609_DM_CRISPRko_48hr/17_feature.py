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

# samples
sample = 'C2'
hue_order = ['GFP', 'mCherry']

df = pd.read_csv('%s/%s/%s_n4_simplified.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
df['r'] = np.sqrt(df['area_nuclear']/math.pi)
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA']/math.pi)
df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt']/df['r']
df_sort = df[df['total_area_ecDNA_sqrt_normalized'] > 0.2].copy().reset_index(drop=True)
df_sample = df_sort[df_sort['group'].isin(hue_order)].copy().reset_index(drop=True)
df = df_sample

df_GFP = df[df['group'] == 'GFP']
df_mCherry = df[df['group'] == 'mCherry']
print(len(df_GFP[df_GFP['n_ecDNA']>10])/len(df_GFP))
print(len(df_mCherry[df_mCherry['n_ecDNA']>10])/len(df_mCherry))


def get_minuslnp(data1: pd.DataFrame, data2: pd.DataFrame, feature: str, limit: int, repeat: int):
    minimum = np.min([len(data1), len(data2)])
    limit = minimum if (minimum < limit) else limit
    p_lst = []
    for j in range(repeat):
        p = -np.log(ks_2samp(data1[feature].sample(n=limit).tolist(), data2[feature].sample(n=limit).tolist())[1])
        p_lst.append(p)
    return limit, np.mean(p_lst)


n_ecDNA_p = get_minuslnp(df_mCherry, df_GFP, 'n_ecDNA', 50, 50)[1]


fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
feature = 'n_ecDNA'
sinaplot(data=df, x='group', y=feature, order=hue_order, violin=False, scale='area', point_size=2)
# plt.ylim([0, 3E8])
plt.savefig('%s/%s/%s_%s_overpoint2.pdf' % (output_dir, sample, sample, feature))
plt.show()