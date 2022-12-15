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

hist_colors = [(0.95, 0.50, 0.50), (0.90, 0.90, 0.90),  (0.50, 0.90, 0.90)]
line_colors = [(0.85, 0.35, 0.25), (0.30, 0.30, 0.30),  (0.30, 0.70, 0.70)]

feature = ['g']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA'])

# g curve
print("Plotting g curve...")
x = np.arange(0, 101, 1)
x_label = 'r'

limit = 80

plt.subplots(figsize=(12, 9))

data = df
number_nuclear = len(data)

mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['g'].tolist())

for i in range(len(data)):
    plt.plot(x[1:limit], data['g'][i][1:limit], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
plt.plot(x[1:limit], mean_curve3[1:limit], color=line_colors[0], label='%s, n=%s' % (sample, number_nuclear))
plt.plot(x[1:limit], ci_lower3[1:limit], color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x[1:limit], ci_higher3[1:limit], color=line_colors[0], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylabel('g')
plt.ylim([0.7, 1.6])
plt.legend()
plt.savefig('%sg_%s.pdf' % (output_dir, sample))
plt.show()

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, x='total_area_ecDNA_sqrt', y='g_value', c=df['MYC_mean'], vmax=40000)
plt.savefig('%s/g_value_vs_total_area_ecDNA_sqrt_%s.pdf' % (output_dir, sample))
plt.show()

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
df['dis_to_hub_area_normalized'] = df['dis_to_hub_area']/np.sqrt(df['area_nuclear'])

xlabel = 'total_area_ecDNA_sqrt'
ylabel = 'g_value'
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