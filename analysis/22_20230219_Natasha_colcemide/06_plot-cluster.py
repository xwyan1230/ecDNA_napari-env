import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
from shared.sinaplot import sinaplot
import math

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230209_analysis_Natasha_DMcolcemid_mchctrl/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 4

# samples
sample1 = 'DMctrl_mchctrl'
figure_name = sample1
pos_threshold = 11000
neg_threshold = 20000
sample1_pos = 'DM H2B-mCherry'
sample1_neg = 'DM'

hue_order = [sample1_neg, sample1_pos]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red

# load data
df1 = pd.read_csv(("%s%s_n%s.txt" % (data_dir, sample1, n_dilation)), na_values=['.'], sep='\t')
df1 = df1[df1['circ_nuclear'] > 0.8].copy().reset_index(drop=True)

feature = ['percentage_area_curve_ecDNA', 'cum_area_ind_ecDNA_filled', 'g']
for f in feature:
    df1[f] = [dat.str_to_float(df1[f][i]) for i in range(len(df1))]

sample_lst = []
for i in range(len(df1)):
    if df1['mCherry_mean'].tolist()[i] < neg_threshold:
        sample_lst.append(sample1_neg)
    elif df1['mCherry_mean'].tolist()[i] > pos_threshold:
        sample_lst.append(sample1_pos)
    else:
        sample_lst.append('NA')
df1['sample'] = sample_lst

df = df1
df['r'] = np.sqrt(df['area_nuclear']/math.pi)
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA']/math.pi)
df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt']/df['r']
df['dis_to_hub_area_normalized'] = df['dis_to_hub_area']/df['r']
df_sample = df[(df['sample'].isin(hue_order)) & (df['circ_nuclear'] > 0.8)].copy().reset_index(drop=True)

# data filter
# df_sort = df_sample
df_sort = df_sample[df_sample['total_area_ecDNA_sqrt_normalized'] > 0.05].copy().reset_index(drop=True)
# df_sort = df_sample[(df_sample['total_area_ecDNA_sqrt_normalized'] > 0.1) & (df_sample['mean_int_DNAFISH'] <= 8500)].copy().reset_index(drop=True)
print(len(df_sort))

# scatter plot
sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_sort, x='total_area_ecDNA_sqrt_normalized', y='dis_to_hub_area_normalized', hue='sample', hue_order=hue_order)
plt.savefig('%s/%s_dis_to_hub_area_normalized_vs_total_area_sqrt_normalized.pdf' % (output_dir, figure_name))
plt.show()

fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
# sns.violinplot(data=df_sort, x='sample', y='dis_to_hub_area_normalized', order=hue_order)
sinaplot(data=df_sort, x='sample', y='dis_to_hub_area_normalized', order=hue_order, violin=False, scale='area')
plt.savefig('%s/%s_dis_to_hub_area_normalized.pdf' % (output_dir, figure_name))
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
plt.savefig('%s%s_dis_barplot.pdf' % (output_dir, figure_name))
plt.show()

# cumulative curve
print("Plotting cumulative curve...")
feature = ['cum_area_ind_ecDNA_filled', 'percentage_area_curve_ecDNA']

for f in feature:
    x_label = 'number of ecDNA hub'
    plt.subplots(figsize=(6, 4))
    for k in range(len(hue_order)):
        data = df_sort[df_sort['sample'] == hue_order[k]].copy().reset_index(drop=True)
        number_nuclear = len(data)

        mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(dat.list_fill_with_last_num(data[f].tolist()))
        for i in range(len(data)):
            plt.plot(data[f][i], alpha=0.05, color=[line_colors[k][j] + 0.05 for j in range(len(line_colors[k]))])
        plt.plot(mean_curve1, color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
        plt.plot(ci_lower1, color=line_colors[k], linestyle='--', linewidth=0.5)
        plt.plot(ci_higher1, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.xlabel(x_label)
    plt.ylabel(f)
    plt.legend()
    plt.savefig('%s/%s_%s.pdf' % (output_dir, figure_name, f))
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
# plt.ylim([0.7, 1.6])
plt.ylim([0.8, 2])
plt.legend()
plt.savefig('%s%s_g.pdf' % (output_dir, figure_name))
plt.show()