import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
from shared.sinaplot import sinaplot
import seaborn as sns
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
sample = 'E9'

GFP_sample = 'GFP'
mCherry_sample = 'mCherry'
hue_order = [GFP_sample, mCherry_sample]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
sns.set_palette(sns.color_palette(line_colors))

df = pd.read_csv('%s/%s/%s_n4_simplified.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
df['total_int_DNAFISH'] = df['mean_int_DNAFISH'] * df['area_nuclear']
feature = ['radial_curve_nuclear', 'radial_curve_DNAFISH', 'percentage_area_curve_ecDNA', 'n_ecDNA_lst',
           'total_area_DNAFISH_lst', 'total_area_ratio_DNAFISH_lst']
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

# radial heatmap
print("Plotting total radial curve...")
column_lst = ['0.025', '0.075', '0.125', '0.175', '0.225', '0.275', '0.325', '0.375', '0.425', '0.475', '0.525',
              '0.575', '0.625', '0.675', '0.725', '0.775', '0.825', '0.875', '0.925', '0.975']
data_heatmap = pd.DataFrame(columns=column_lst)
for s in hue_order:
    data_sample = df[df['group'] == s].copy().reset_index(drop=True)
    data_radial = pd.DataFrame()
    for j in range(len(column_lst)):
        data_radial[column_lst[j]] = [data_sample['radial_curve_DNAFISH'][i][j] for i in range(len(data_sample))]
    data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()
data_heatmap.index = hue_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/%s/%s_heatmap_DNAFISH.pdf' % (output_dir, sample, sample))
plt.close()

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=1.3, vmin=0.6, square=True, cmap='coolwarm')
plt.savefig('%s/%s/%s_heatmap_DNAFISH_fixscale.pdf' % (output_dir, sample, sample))
plt.close()

# radial curve
x = np.arange(0.025, 1, 0.05)
x_label = 'relative r'
plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['group'] == hue_order[k]].copy().reset_index(drop=True)
    number_nuclear = len(data)
    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_DNAFISH'].tolist())
    for i in range(len(data)):
        plt.plot(x, data['radial_curve_DNAFISH'][i], alpha=0.01, color=line_colors[k])
    plt.plot(x, mean_curve3, color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[k], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.4, 1.6])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s/%s_radial_curve_DNAFISH.pdf' % (output_dir, sample, sample))
plt.close()

# n_ecDNA
print("Plotting single feature...")
fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
feature = 'n_ecDNA'
sinaplot(data=df, x='group', y=feature, order=hue_order, violin=False, scale='area', point_size=2)
plt.ylim([0, 20])
plt.savefig('%s/%s/%s_%s.pdf' % (output_dir, sample, sample, feature))
plt.close()

fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
feature = 'total_area_ecDNA'
sinaplot(data=df, x='group', y=feature, order=hue_order, violin=False, scale='area', point_size=2)
plt.ylim([0, 5000])
plt.savefig('%s/%s/%s_%s.pdf' % (output_dir, sample, sample, feature))
plt.close()

fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
feature = 'total_area_ratio_ecDNA'
sinaplot(data=df, x='group', y=feature, order=hue_order, violin=False, scale='area', point_size=2)
plt.ylim([0, 0.3])
plt.savefig('%s/%s/%s_%s.pdf' % (output_dir, sample, sample, feature))
plt.close()

fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
feature = 'mean_int_DNAFISH'
sinaplot(data=df, x='group', y=feature, order=hue_order, violin=False, scale='area', point_size=2)
plt.ylim([0, 15000])
plt.savefig('%s/%s/%s_%s.pdf' % (output_dir, sample, sample, feature))
plt.close()

fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
feature = 'total_int_DNAFISH'
sinaplot(data=df, x='group', y=feature, order=hue_order, violin=False, scale='area', point_size=2)
plt.ylim([0, 3E8])
plt.savefig('%s/%s/%s_%s.pdf' % (output_dir, sample, sample, feature))
plt.close()

# CPF
print("Plotting cpf...")
plt.subplots(figsize=(12, 9))
feature = 'n_ecDNA'
histrange = [0, 20]
plt.hist(df_GFP[feature].tolist(), 40, density=True, histtype='step', cumulative=True, label='GFP, n=%s' % len(df_GFP), range=histrange, color=line_colors[0])
plt.hist(df_mCherry[feature].tolist(), 40, density=True, histtype='step', cumulative=True, label='mCherry, n=%s' % len(df_mCherry), range=histrange, color=line_colors[1])
plt.xlabel(feature)
plt.ylabel('CPF')  # cumulative distribution function
plt.legend()
plt.savefig('%s/%s/%s_cpf_%s.pdf' % (output_dir, sample, sample, feature))
plt.close()

plt.subplots(figsize=(12, 9))
feature = 'total_area_ratio_ecDNA'
histrange = [0, 0.2]
plt.hist(df_GFP[feature].tolist(), 40, density=True, histtype='step', cumulative=True, label='GFP, n=%s' % len(df_GFP), range=histrange, color=line_colors[0])
plt.hist(df_mCherry[feature].tolist(), 40, density=True, histtype='step', cumulative=True, label='mCherry, n=%s' % len(df_mCherry), range=histrange, color=line_colors[1])
plt.xlabel(feature)
plt.ylabel('CPF')  # cumulative distribution function
plt.legend()
plt.savefig('%s/%s/%s_cpf_%s.pdf' % (output_dir, sample, sample, feature))
plt.close()

plt.subplots(figsize=(12, 9))
feature = 'total_area_ecDNA'
histrange = [0, 2500]
plt.hist(df_GFP[feature].tolist(), 40, density=True, histtype='step', cumulative=True, label='GFP, n=%s' % len(df_GFP), range=histrange, color=line_colors[0])
plt.hist(df_mCherry[feature].tolist(), 40, density=True, histtype='step', cumulative=True, label='mCherry, n=%s' % len(df_mCherry), range=histrange, color=line_colors[1])
plt.xlabel(feature)
plt.ylabel('CPF')  # cumulative distribution function
plt.legend()
plt.savefig('%s/%s/%s_cpf_%s.pdf' % (output_dir, sample, sample, feature))
plt.close()

plt.subplots(figsize=(12, 9))
feature = 'mean_int_DNAFISH'
histrange = [0, 15000]
plt.hist(df_GFP[feature].tolist(), 40, density=True, histtype='step', cumulative=True, label='GFP, n=%s' % len(df_GFP), range=histrange, color=line_colors[0])
plt.hist(df_mCherry[feature].tolist(), 40, density=True, histtype='step', cumulative=True, label='mCherry, n=%s' % len(df_mCherry), range=histrange, color=line_colors[1])
plt.xlabel(feature)
plt.ylabel('CPF')  # cumulative distribution function
plt.legend()
plt.savefig('%s/%s/%s_cpf_%s.pdf' % (output_dir, sample, sample, feature))
plt.close()

plt.subplots(figsize=(12, 9))
feature = 'total_int_DNAFISH'
histrange = [0, 3.5E8]
plt.hist(df_GFP[feature].tolist(), 40, density=True, histtype='step', cumulative=True, label='GFP, n=%s' % len(df_GFP), range=histrange, color=line_colors[0])
plt.hist(df_mCherry[feature].tolist(), 40, density=True, histtype='step', cumulative=True, label='mCherry, n=%s' % len(df_mCherry), range=histrange, color=line_colors[1])
plt.xlabel(feature)
plt.ylabel('CPF')  # cumulative distribution function
plt.legend()
plt.savefig('%s/%s/%s_cpf_%s.pdf' % (output_dir, sample, sample, feature))
plt.close()

# cum_curve
print("Plotting total cum curve...")
f = 'percentage_area_curve_ecDNA'
x_label = 'number of ecDNA hub'
df1 = df_GFP
df2 = df_mCherry
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
plt.savefig('%s/%s/%s_%s.pdf' % (output_dir, sample, sample, f))
plt.close()

plt.subplots(figsize=(6, 4))
for i in range(len(df1)):
    plt.plot(df1[f][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
plt.plot(mean_curve1, color=line_colors[0], label='%s, n=%s' % (hue_order[0], number_nuclear1))
plt.plot(ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylabel(f)
plt.legend()
plt.savefig('%s/%s/%s_%s_%s.pdf' % (output_dir, sample, sample, f, hue_order[0]))
plt.close()

plt.subplots(figsize=(6, 4))
for i in range(len(df2)):
    plt.plot(df2[f][i], alpha=0.05, color=[line_colors[1][j]+0.05 for j in range(len(line_colors[1]))])
plt.plot(mean_curve2, color=line_colors[1], label='%s, n=%s' % (hue_order[1], number_nuclear2))
plt.plot(ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylabel(f)
plt.legend()
plt.savefig('%s/%s/%s_%s_%s.pdf' % (output_dir, sample, sample, f, hue_order[1]))
plt.close()

# n_ecDNA, intensity distribution curve
print("Plotting n_ecDNA, intensity distribution...")
x = np.arange(2500, 15*5000+2500, 5000)
x_label = 'intensity'
feature = 'n_ecDNA_lst'
plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['group'] == hue_order[k]].copy().reset_index(drop=True)
    number_nuclear = len(data)
    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data[feature].tolist())
    for i in range(len(data)):
        plt.plot(x, data[feature][i], alpha=0.001, color=line_colors[k])
    plt.plot(x, mean_curve3, color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[k], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylim([0, 20])
plt.ylabel(feature)
plt.legend()
plt.savefig('%s%s/%s_%s.pdf' % (output_dir, sample, sample, feature))
plt.close()

x = np.arange(2500, 15*5000+2500, 5000)
x_label = 'intensity'
feature = 'total_area_DNAFISH_lst'
plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['group'] == hue_order[k]].copy().reset_index(drop=True)
    number_nuclear = len(data)
    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data[feature].tolist())
    for i in range(len(data)):
        plt.plot(x, data[feature][i], alpha=0.005, color=line_colors[k])
    plt.plot(x, mean_curve3, color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[k], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylim([0, 20000])
plt.ylabel(feature)
plt.legend()
plt.savefig('%s%s/%s_%s.pdf' % (output_dir, sample, sample, feature))
plt.close()

x = np.arange(2500, 15*5000+2500, 5000)
x_label = 'intensity'
feature = 'total_area_ratio_DNAFISH_lst'
plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['group'] == hue_order[k]].copy().reset_index(drop=True)
    number_nuclear = len(data)
    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data[feature].tolist())
    for i in range(len(data)):
        plt.plot(x, data[feature][i], alpha=0.005, color=line_colors[k])
    plt.plot(x, mean_curve3, color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[k], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylim([0, 1])
plt.ylabel(feature)
plt.legend()
plt.savefig('%s%s/%s_%s.pdf' % (output_dir, sample, sample, feature))
plt.close()

### SORT
df_sort = df[df['total_area_ecDNA_sqrt_normalized'] > 0.2].copy().reset_index(drop=True)
df = df_sort
df_GFP = df[df['group'] == 'GFP'].copy().reset_index(drop=True)
df_mCherry = df[df['group'] == 'mCherry'].copy().reset_index(drop=True)

# radial heatmap
print("Plotting sort radial curve...")
column_lst = ['0.025', '0.075', '0.125', '0.175', '0.225', '0.275', '0.325', '0.375', '0.425', '0.475', '0.525',
              '0.575', '0.625', '0.675', '0.725', '0.775', '0.825', '0.875', '0.925', '0.975']
data_heatmap = pd.DataFrame(columns=column_lst)
for s in hue_order:
    data_sample = df[df['group'] == s].copy().reset_index(drop=True)
    data_radial = pd.DataFrame()
    for j in range(len(column_lst)):
        data_radial[column_lst[j]] = [data_sample['radial_curve_DNAFISH'][i][j] for i in range(len(data_sample))]
    data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()
data_heatmap.index = hue_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/%s/%s_heatmap_DNAFISH_point2.pdf' % (output_dir, sample, sample))
plt.close()

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=1.6, vmin=0.4, square=True, cmap='coolwarm')
plt.savefig('%s/%s/%s_heatmap_DNAFISH_fixscale_point2.pdf' % (output_dir, sample, sample))
plt.close()

# radial curve
x = np.arange(0.025, 1, 0.05)
x_label = 'relative r'
plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['group'] == hue_order[k]].copy().reset_index(drop=True)
    number_nuclear = len(data)
    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_DNAFISH'].tolist())
    for i in range(len(data)):
        plt.plot(x, data['radial_curve_DNAFISH'][i], alpha=0.01, color=line_colors[k])
    plt.plot(x, mean_curve3, color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[k], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.4, 1.6])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s/%s_radial_curve_DNAFISH_point2.pdf' % (output_dir, sample, sample))
plt.close()

# n_ecDNA
print("Plotting sort n_ecDNA...")
fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
feature = 'n_ecDNA'
sinaplot(data=df, x='group', y=feature, order=hue_order, violin=False, scale='area', point_size=2)
plt.ylim([0, 20])
plt.savefig('%s/%s/%s_%s_point2.pdf' % (output_dir, sample, sample, feature))
plt.close()

# averageD
print("Plotting sort averageD...")
fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
sinaplot(data=df, x='group', y='dis_to_hub_area_normalized', order=hue_order, violin=False, scale='area')
plt.savefig('%s/%s/%s_dis_to_hub_area_normalized_point2.pdf' % (output_dir, sample, sample))
plt.close()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_GFP, y='dis_to_hub_area_normalized', x='total_area_ecDNA_sqrt_normalized', label='GFP, n=%s' % len(df_GFP), color=line_colors[0], s=10, alpha=1)
sns.scatterplot(data=df_mCherry, y='dis_to_hub_area_normalized', x='total_area_ecDNA_sqrt_normalized', label='mCherry, n=%s' % len(df_mCherry), color=line_colors[1], s=10, alpha=1)
plt.legend()
plt.savefig('%s/%s/%s_dis_to_hub_area_normalized_scatter_point2.pdf' % (output_dir, sample, sample))
plt.close()

# cum_curve
print("Plotting sort cum curve...")
f = 'percentage_area_curve_ecDNA'
x_label = 'number of ecDNA hub'
df1 = df_GFP
df2 = df_mCherry
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
plt.savefig('%s/%s/%s_%s_point2.pdf' % (output_dir, sample, sample, f))
plt.close()

plt.subplots(figsize=(6, 4))
for i in range(len(df1)):
    plt.plot(df1[f][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
plt.plot(mean_curve1, color=line_colors[0], label='%s, n=%s' % (hue_order[0], number_nuclear1))
plt.plot(ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylabel(f)
plt.legend()
plt.savefig('%s/%s/%s_%s_%s_point2.pdf' % (output_dir, sample, sample, f, hue_order[0]))
plt.close()

plt.subplots(figsize=(6, 4))
for i in range(len(df2)):
    plt.plot(df2[f][i], alpha=0.05, color=[line_colors[1][j]+0.05 for j in range(len(line_colors[1]))])
plt.plot(mean_curve2, color=line_colors[1], label='%s, n=%s' % (hue_order[1], number_nuclear2))
plt.plot(ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylabel(f)
plt.legend()
plt.savefig('%s/%s/%s_%s_%s_point2.pdf' % (output_dir, sample, sample, f, hue_order[1]))
plt.close()

print("DONE!")