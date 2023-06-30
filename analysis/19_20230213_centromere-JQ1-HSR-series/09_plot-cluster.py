import skimage.io as skio
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import math
from sklearn.linear_model import LinearRegression
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230206_analysis_centromere-JQ1-DM-HSR/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'DM-mCh_mix_HSR'
pos_threshold = 10000
neg_threshold = 5500
pos = 'DM H2B-mCherry'
neg = 'HSR'
hue_order = [pos, neg]
n_dilation = 4

df = pd.read_csv(("%s%s_n%s.txt" % (data_dir, sample, n_dilation)), na_values=['.'], sep='\t')

hist_colors = [(0.90, 0.90, 0.90), (0.95, 0.50, 0.50),  (0.50, 0.90, 0.90)]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25),  (0.30, 0.70, 0.70)]
# '#FFA500', '#40E0D0', '#DA70D6', '#00CED1'
feature = ['percentage_area_curve_ecDNA', 'cum_area_ind_ecDNA_filled', 'g']

for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

sample_lst = []
for i in range(len(df)):
    if df['mCherry_mean'].tolist()[i] < neg_threshold:
        sample_lst.append(neg)
    elif df['mCherry_mean'].tolist()[i] > pos_threshold:
        sample_lst.append(pos)
    else:
        sample_lst.append('NA')
df['sample'] = sample_lst
df['r'] = np.sqrt(df['area_nuclear']/math.pi)
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA']/math.pi)
df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt']/df['r']
df['dis_to_hub_area_normalized'] = df['dis_to_hub_area']/df['r']

# df_sort = df[df['sample'].isin([neg, pos])].copy().reset_index()
# df_sort = df[(df['sample'].isin([neg, pos])) & (df['circ_nuclear'] > 0.8)].copy().reset_index()
# df_sort = df[(df['sample'].isin([neg, pos])) & (df['dis_to_hub_area_normalized'] < 1.2)].copy().reset_index(drop=True)
df_sort = df[(df['sample'].isin([neg, pos])) & (df['dis_to_hub_area_normalized'] < 1.2) & (df['total_area_ecDNA_sqrt_normalized'] > 0.2)].copy().reset_index(drop=True)
df_pos = df_sort[df_sort['sample'].isin([pos])].copy().reset_index(drop=True)
df_neg = df_sort[df_sort['sample'].isin([neg])].copy().reset_index(drop=True)
df_sort['total_area_ecDNA_sqrt'] = np.sqrt(df_sort['total_area_ecDNA'])
seq = [df_pos, df_neg]
print(len(df_sort))

# scatter plot
sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_sort, x='total_area_ecDNA_sqrt_normalized', y='dis_to_hub_area_normalized', hue='sample', hue_order=hue_order)
plt.savefig('%s/%s_dis_to_hub_area_normalized_vs_total_area_sqrt_normalized.pdf' % (output_dir, sample))
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_sort, x='total_area_ecDNA_sqrt_normalized', y='dis_to_hub_area_normalized', c=df_sort['circ_nuclear'].tolist())
plt.savefig('%s/%s_dis_to_hub_area_normalized_vs_total_area_sqrt_normalized_color_by_circ.pdf' % (output_dir, sample))
plt.show()

fig, ax = plt.subplots(figsize=(5, 9))
fig.subplots_adjust(left=0.2)
sns.violinplot(data=df_sort, x='sample', y='dis_to_hub_area_normalized', order=hue_order)
plt.savefig('%s/%s_dis_to_hub_area_normalized.pdf' % (output_dir, sample))
plt.show()

d_range = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
mean_df = []
for i in range(len(hue_order)):
    data_sample = df_sort[df_sort['sample'] == hue_order[i]].copy().reset_index(drop=True)
    temp = []
    for j in range(len(d_range)-1):
        temp.append(len(data_sample[(data_sample['dis_to_hub_area_normalized'] >= d_range[j]) & (data_sample['dis_to_hub_area_normalized'] < d_range[j+1])]) * 1.0 / len(data_sample))
    mean_df.append(temp)

df_temp = pd.DataFrame(mean_df)
df_temp.columns = ['%.1f' % elem for elem in d_range][1:]
df_temp.index = hue_order

dis = pd.DataFrame()
dis['sample'] = list(hue_order) * 6
dis['dis_to_hub'] = [0.1] * len(hue_order) + [0.3] * len(hue_order) + [0.5] * len(hue_order) + [0.7] * len(hue_order) + [0.9] * len(hue_order) + [1.1] * len(hue_order)
dis['percentage'] = df_temp['0.2'].tolist() + df_temp['0.4'].tolist() + df_temp['0.6'].tolist() + df_temp['0.8'].tolist() + df_temp['1.0'].tolist() + df_temp['1.2'].tolist()
print(dis.head())

plt.subplots(figsize=(6, 9))
sns.barplot(data=dis, x='sample', y='percentage', hue='dis_to_hub')
plt.savefig('%s%s_dis_barplot.pdf' % (output_dir, sample))
plt.show()

# cumulative curve
print("Plotting cumulative curve...")
feature = ['cum_area_ind_ecDNA_filled', 'percentage_area_curve_ecDNA']

for f in feature:
    x_label = 'number of ecDNA hub'

    df1 = seq[0]
    df2 = seq[1]
    number_nuclear1 = len(df1)
    number_nuclear2 = len(df2)

    mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(dat.list_fill_with_last_num(df1[f].tolist()))
    mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(dat.list_fill_with_last_num(df2[f].tolist()))

    plt.subplots(figsize=(6, 4))
    for i in range(len(df1)):
        plt.plot(df1[f][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
    for i in range(len(df2)):
        plt.plot(df2[f][i], alpha=0.05, color=[line_colors[1][j]+0.05 for j in range(len(line_colors[1]))])
    plt.plot(mean_curve1, color=line_colors[0], label='%s, n=%s' % (hue_order[0], number_nuclear1))
    plt.plot(mean_curve2, color=line_colors[1], label='%s, n=%s' % (hue_order[1], number_nuclear2))
    plt.plot(ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.xlabel(x_label)
    plt.ylabel(f)
    plt.legend()
    plt.savefig('%s/%s_%s.pdf' % (output_dir, sample, f))
    plt.show()

# g
print("Plotting g curve...")
x = np.arange(0, 101, 1)
x_label = 'r'

limit = 80

plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df_sort[df_sort['sample'] == hue_order[k]].copy().reset_index(drop=True)
    number_nuclear = len(data)

    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['g'].tolist())

    for i in range(len(data)):
        plt.plot(x[1:limit], data['g'][i][1:limit], alpha=0.05, color=[line_colors[k][j] + 0.05 for j in range(len(line_colors[k]))])
    plt.plot(x[1:limit], mean_curve3[1:limit], color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x[1:limit], ci_lower3[1:limit], color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x[1:limit], ci_higher3[1:limit], color=line_colors[k], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylabel('g')
plt.ylim([0.7, 1.6])
# plt.ylim([0.8, 2])
plt.legend()
plt.savefig('%s%s_g.pdf' % (output_dir, sample))
plt.show()