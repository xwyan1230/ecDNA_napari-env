import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import matplotlib as mpl
import numpy as np
from matplotlib import cm
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221215_analysis_MYC-IF/20221117_immunoFISH_acid/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'MYC-old'
n = 8
df = pd.read_csv(("%s%s_radial_n%s.txt" % (data_dir, sample, n)), na_values=['.'], sep='\t')

hist_colors = [(0.90, 0.90, 0.90), (0.95, 0.50, 0.50),  (0.50, 0.90, 0.90)]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25),  (0.30, 0.70, 0.70)]

feature = ['radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_normalized']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

# radial curve
print("Plotting radial curve...")
x = np.arange(0.05, 1, 0.1)
x_label = 'relative r'

sample_lst = []
for i in range(len(df)):
    if df['MYC_mean'].tolist()[i] < 10000:
        sample_lst.append(10000)
    elif df['MYC_mean'].tolist()[i] < 15000:
        sample_lst.append(15000)
    elif df['MYC_mean'].tolist()[i] < 20000:
        sample_lst.append(20000)
    elif df['MYC_mean'].tolist()[i] < 25000:
        sample_lst.append(25000)
    elif df['MYC_mean'].tolist()[i] < 30000:
        sample_lst.append(30000)
    elif df['MYC_mean'].tolist()[i] < 35000:
        sample_lst.append(35000)
    elif df['MYC_mean'].tolist()[i] < 40000:
        sample_lst.append(40000)
    elif df['MYC_mean'].tolist()[i] < 45000:
        sample_lst.append(45000)
    else:
        sample_lst.append(50000)

df['sample'] = sample_lst

hue_order = [10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000]

df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA'])

xlabel = 'total_area_ecDNA_sqrt'
ylabel = 'MYC_mean'
plt.subplots(figsize=(12, 9))
norm = mpl.colors.Normalize(vmin=10000, vmax=65000)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.Spectral_r)
for i in range(len(hue_order)):
    data = df[df['sample'] == hue_order[i]].copy().reset_index(drop=True)
    if len(data) > 0:
        number_nuclear = len(data)
        plt.scatter(data[xlabel], data[ylabel], color=mapper.to_rgba(hue_order)[i], label='%s, n=%s' % (hue_order[i], number_nuclear))
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend()
plt.savefig('%s%s_vs_%s_%s_by_MYC_group.pdf' % (output_dir, ylabel, xlabel, sample))
plt.show()

