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

sample = 'MYC-new'
df = pd.read_csv(("%s%s_ecDNA.txt" % (data_dir, sample)), na_values=['.'], sep='\t')

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

plt.subplots(figsize=(12, 9))
norm = mpl.colors.Normalize(vmin=10000, vmax=65000)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.Spectral_r)
for i in range(len(hue_order)):
    data = df[df['sample'] == hue_order[i]].copy().reset_index(drop=True)
    if len(data) > 0:
        number_nuclear = len(data)
        mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_normalized'].tolist())
        plt.plot(x, mean_curve3, color=mapper.to_rgba(hue_order)[i], label='%s, n=%s' % (hue_order[i], number_nuclear))
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.7, 1.3])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%sradial_curve_normalized_by_MYC_group_%s.pdf' % (output_dir, sample))
plt.show()

data_heatmap = pd.DataFrame(columns=['0.05', '0.15', '0.25', '0.35', '0.45', '0.55', '0.65', '0.75', '0.85', '0.95'])

for s in hue_order:
    data_sample = df[df['sample'] == s].copy().reset_index(drop=True)
    data_radial = pd.DataFrame()
    data_radial['0.05'] = [data_sample['radial_curve_normalized'][i][0] for i in range(len(data_sample))]
    data_radial['0.15'] = [data_sample['radial_curve_normalized'][i][1] for i in range(len(data_sample))]
    data_radial['0.25'] = [data_sample['radial_curve_normalized'][i][2] for i in range(len(data_sample))]
    data_radial['0.35'] = [data_sample['radial_curve_normalized'][i][3] for i in range(len(data_sample))]
    data_radial['0.45'] = [data_sample['radial_curve_normalized'][i][4] for i in range(len(data_sample))]
    data_radial['0.55'] = [data_sample['radial_curve_normalized'][i][5] for i in range(len(data_sample))]
    data_radial['0.65'] = [data_sample['radial_curve_normalized'][i][6] for i in range(len(data_sample))]
    data_radial['0.75'] = [data_sample['radial_curve_normalized'][i][7] for i in range(len(data_sample))]
    data_radial['0.85'] = [data_sample['radial_curve_normalized'][i][8] for i in range(len(data_sample))]
    data_radial['0.95'] = [data_sample['radial_curve_normalized'][i][9] for i in range(len(data_sample))]

    data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()

data_heatmap.index = hue_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_%s.pdf' % (output_dir, sample))
plt.show()