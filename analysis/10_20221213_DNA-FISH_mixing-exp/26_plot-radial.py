import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221213_analysis_DNA-FISH_mixing_exp/20221121_H2B-series_POLR3D-BRD4-BRD1-DAPK2/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 4

# samples
sample1 = 'H2B+POLR3D'
figure_name = sample1
pos_threshold = 20000
neg_threshold = 12000
sample1_pos = 'DM H2B-mCherry'
sample1_neg = 'DM POLR3Dko'

hue_order = [sample1_pos, sample1_neg]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
hue_order1 = ['background', 'DNAFISH']

# load data
df1 = pd.read_csv(("%s%s_n%s.txt" % (data_dir, sample1, n_dilation)), na_values=['.'], sep='\t')

feature = ['radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_normalized', 'radial_curve_DNAFISH_seg',
           'nuclear_int', 'DNAFISH_int', 'DNAFISH_seg_label', 'int_r_to_edge', 'int_r_to_center', 'int_relative_r']
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
df_sample = df[df['sample'].isin(hue_order)].copy().reset_index(drop=True)
df = df_sample


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask) if e != 0]
    out = list(map(img.__getitem__, index))
    return out


# heatmap
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
plt.savefig('%s/heatmap_normalized_n%s_%s.pdf' % (output_dir, n_dilation, figure_name))
plt.show()

data_heatmap = pd.DataFrame(columns=['0.05', '0.15', '0.25', '0.35', '0.45', '0.55', '0.65', '0.75', '0.85', '0.95'])

for s in hue_order:
    data_sample = df[df['sample'] == s].copy().reset_index(drop=True)
    data_radial = pd.DataFrame()
    data_radial['0.05'] = [data_sample['radial_curve_DNAFISH'][i][0] for i in range(len(data_sample))]
    data_radial['0.15'] = [data_sample['radial_curve_DNAFISH'][i][1] for i in range(len(data_sample))]
    data_radial['0.25'] = [data_sample['radial_curve_DNAFISH'][i][2] for i in range(len(data_sample))]
    data_radial['0.35'] = [data_sample['radial_curve_DNAFISH'][i][3] for i in range(len(data_sample))]
    data_radial['0.45'] = [data_sample['radial_curve_DNAFISH'][i][4] for i in range(len(data_sample))]
    data_radial['0.55'] = [data_sample['radial_curve_DNAFISH'][i][5] for i in range(len(data_sample))]
    data_radial['0.65'] = [data_sample['radial_curve_DNAFISH'][i][6] for i in range(len(data_sample))]
    data_radial['0.75'] = [data_sample['radial_curve_DNAFISH'][i][7] for i in range(len(data_sample))]
    data_radial['0.85'] = [data_sample['radial_curve_DNAFISH'][i][8] for i in range(len(data_sample))]
    data_radial['0.95'] = [data_sample['radial_curve_DNAFISH'][i][9] for i in range(len(data_sample))]

    data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()

data_heatmap.index = hue_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_DNAFISH_n%s_%s.pdf' % (output_dir, n_dilation, figure_name))
plt.show()

# radial curve
print("Plotting radial curve...")
x = np.arange(0.05, 1, 0.1)
x_label = 'relative r'

for s in hue_order:
    data = df[df['sample'] == s].copy().reset_index(drop=True)
    number_nuclear = len(data)

    mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data['radial_curve_DNAFISH'].tolist())
    mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(data['radial_curve_nuclear'].tolist())
    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_normalized'].tolist())

    plt.subplots(figsize=(12, 9))
    for i in range(len(data)):
        plt.plot(x, data['radial_curve_DNAFISH'][i], alpha=0.05, color=[line_colors[1][j] + 0.05 for j in range(len(line_colors[1]))])
    for i in range(len(data)):
        plt.plot(x, data['radial_curve_nuclear'][i], alpha=0.05, color=[line_colors[0][j]+0.05 for j in range(len(line_colors[0]))])
    plt.plot(x, mean_curve1, color=line_colors[1], label='%s, n=%s' % ('ecDNA (MYC DNA FISH)', number_nuclear))
    plt.plot(x, mean_curve2, color=line_colors[0], label='%s, n=%s' % ('DNA (hoechst stain)', number_nuclear))
    plt.plot(x, ci_lower1, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher1, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_lower2, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher2, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.xlabel(x_label)
    plt.ylim([0.5, 1.5])
    plt.ylabel('radial_curve')
    plt.legend()
    plt.savefig('%sradial_curve_n%s_%s_%s.pdf' % (output_dir, n_dilation, figure_name, s))
    plt.show()

"""    plt.subplots(figsize=(12, 9))
    for i in range(len(data)):
        plt.plot(x, data['radial_curve_normalized'][i], alpha=0.05, color=[line_colors[1][j] + 0.05 for j in range(len(line_colors[1]))])
    plt.plot(x, mean_curve3, color=line_colors[1], label='%s, n=%s' % ('ecDNA (normalized)', number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.axhline(y=1, color='black', linestyle='--')
    plt.xlabel(x_label)
    plt.ylim([0.5, 1.5])
    plt.ylabel('radial_curve')
    plt.legend()
    plt.savefig('%sradial_curve_normalized_n%s_%s_%s.pdf' % (output_dir, n_dilation, figure_name, s))
    plt.show()"""

plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['sample'] == hue_order[k]].copy().reset_index(drop=True)
    number_nuclear = len(data)

    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_normalized'].tolist())

    for i in range(len(data)):
        plt.plot(x, data['radial_curve_normalized'][i], alpha=0.05, color=[line_colors[k][j] + 0.05 for j in range(len(line_colors[k]))])
    plt.plot(x, mean_curve3, color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[k], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.5, 1.5])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%sradial_curve_normalized_n%s_%s.pdf' % (output_dir, n_dilation, figure_name))
plt.show()

plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['sample'] == hue_order[k]].copy().reset_index(drop=True)
    number_nuclear = len(data)

    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_DNAFISH'].tolist())

    for i in range(len(data)):
        plt.plot(x, data['radial_curve_DNAFISH'][i], alpha=0.05, color=[line_colors[k][j] + 0.05 for j in range(len(line_colors[k]))])
    plt.plot(x, mean_curve3, color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[k], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.5, 1.5])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%sradial_curve_DNAFISH_n%s_%s.pdf' % (output_dir, n_dilation, figure_name))
plt.show()

# new approach
df_r = pd.DataFrame()
for k in range(len(hue_order)):
    data = df[df['sample'] == hue_order[k]].copy().reset_index(drop=True)
    df_r_temp = pd.DataFrame()
    sample_category_lst = []
    intensity_lst = []
    seg_lst = []
    relative_r_lst = []
    for i in range(len(data)):
        sample_category_lst = sample_category_lst + [hue_order1[0]] * len(data['int_relative_r'][i]) + [hue_order1[1]] * len(data['int_relative_r'][i])
        relative_r_lst = relative_r_lst + data['int_relative_r'][i] + data['int_relative_r'][i]
        intensity_lst = intensity_lst + data['nuclear_int'][i] + data['DNAFISH_int'][i]
        seg_lst = seg_lst + [1] * len(data['int_relative_r'][i]) + data['DNAFISH_seg_label'][i]
    df_r_temp['sample_category'] = sample_category_lst
    df_r_temp['relative_r'] = relative_r_lst
    df_r_temp['intensity'] = intensity_lst
    df_r_temp['seg'] = seg_lst
    df_r_temp['sample'] = [hue_order[k]] * len(df_r_temp)
    len_bg = len(df_r_temp[df_r_temp['sample_category'] == hue_order1[0]])
    len_sample = len(df_r_temp[(df_r_temp['sample_category'] == hue_order1[1]) & (df_r_temp['seg'] == 1)])
    df_r_temp['weights'] = [1.0/len_bg if i == hue_order1[0] else 1.0/len_sample for i in df_r_temp['sample_category']]
    df_r = pd.concat([df_r, df_r_temp], axis=0)

for k in range(len(hue_order)):
    data = df[df['sample'] == hue_order[k]].copy().reset_index(drop=True)
    data_r = df_r[df_r['sample'] == hue_order[k]].copy().reset_index(drop=True)
    sns.set_palette(sns.color_palette(line_colors))

    data_r_sort = data_r[data_r['seg'] == 1].copy().reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(9, 6))
    fig.subplots_adjust(left=0.2)
    ax = sns.histplot(data=data_r_sort, x='relative_r', hue='sample_category', hue_order=hue_order1, multiple='dodge', bins=20, weights=data_r_sort['weights'])
    plt.savefig('%s/%s_histplot_DNAFISH_seg_n%s_cell%s_%s.pdf' % (output_dir, figure_name, n_dilation, len(data), hue_order[k]))
    plt.show()

# heatmap
x = np.arange(0.025, 1, 0.05)
x_label = 'relative r'
up = 20
df_seg_normalized = pd.DataFrame(columns=x[:up])
df_int_normalized = pd.DataFrame(columns=x[:up])
df_int_seg_normalized = pd.DataFrame(columns=x[:up])

# seg
for k in range(len(hue_order)):
    data = df[df['sample'] == hue_order[k]].copy().reset_index(drop=True)
    data_r = df_r[df_r['sample'] == hue_order[k]].copy().reset_index(drop=True)

    data_heatmap = pd.DataFrame(columns=x[:up])
    for s in hue_order1:
        data_sample = data_r[(data_r['sample_category'] == s) & (data_r['seg'] == 1)].copy().reset_index(drop=True)
        total_data_sample = len(data_sample)
        data_radial = []
        for i in range(len(x[:up])):
            if x[i] == 0.025:
                n_data_radial = len(data_sample[(data_sample['relative_r'] >= 0) & (data_sample['relative_r'] <= 0.05)])
            else:
                n_data_radial = len(data_sample[(data_sample['relative_r'] > (x[i]-0.025)) & (data_sample['relative_r'] <= (x[i]+0.025))])
            data_radial.append(n_data_radial * 1.0/total_data_sample)
        data_heatmap.loc[len(data_heatmap.index)] = data_radial
    data_heatmap.index = hue_order1
    data_heatmap.columns = ['%.3f' % elem for elem in x[:up]]

    plt.subplots(figsize=(12, len(hue_order1)))
    ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
    plt.savefig('%s/%s_heatmap_DNAFISH_seg_before_normalization_n%s_cell%s_%s.pdf' % (output_dir, figure_name, n_dilation, len(data), hue_order[k]))
    plt.show()

    data_heatmap_normalized = pd.DataFrame(columns=x[:up])
    temp = []
    for i in data_heatmap.columns:
        temp.append(data_heatmap[i][hue_order1[1]] / data_heatmap[i][hue_order1[0]])
    data_heatmap_normalized.loc[0] = temp
    data_heatmap_normalized.index = [hue_order1[1]]
    data_heatmap_normalized.columns = ['%.3f' % elem for elem in x[:up]]
    df_seg_normalized.loc[len(df_seg_normalized.index)] = temp

    plt.subplots(figsize=(12, len(hue_order)))
    ax1 = sns.heatmap(data_heatmap_normalized, cbar=0, linewidths=2, vmax=data_heatmap_normalized.values.max(),
                      vmin=data_heatmap_normalized.values.min(), square=True, cmap='coolwarm')
    plt.savefig('%s/%s_heatmap_DNAFISH_seg_n%s_cell%s_%s.pdf' % (output_dir, figure_name, n_dilation, len(data), hue_order[k]))
    plt.show()

# int
for k in range(len(hue_order)):
    data = df[df['sample'] == hue_order[k]].copy().reset_index(drop=True)
    data_r = df_r[df_r['sample'] == hue_order[k]].copy().reset_index(drop=True)
    data_heatmap = pd.DataFrame(columns=x[:up])

    for s in hue_order1:
        data_sample = data_r[data_r['sample_category'] == s].copy().reset_index(drop=True)
        data_radial = []
        for i in range(len(x[:up])):
            if x[i] == 0.025:
                filter = (data_sample['relative_r'] >= 0) & (data_sample['relative_r'] <= 0.05)
            else:
                filter = (data_sample['relative_r'] > (x[i]-0.025)) & (data_sample['relative_r'] <= (x[i]+0.025))
            data_radial.append(sum(data_sample[filter]['intensity'])/sum(data_sample['intensity']))
        data_heatmap.loc[len(data_heatmap.index)] = data_radial
    data_heatmap.index = hue_order1
    data_heatmap.columns = ['%.3f' % elem for elem in x[:up]]

    plt.subplots(figsize=(12, len(hue_order1)))
    ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
    plt.savefig('%s/%s_heatmap_DNAFISH_int_before_normalization_n%s_cell%s_%s.pdf' % (output_dir, figure_name, n_dilation, len(data), hue_order[k]))
    plt.show()

    data_heatmap_normalized = pd.DataFrame(columns=x[:up])
    temp = []
    for i in data_heatmap.columns:
        temp.append(data_heatmap[i][hue_order1[1]] / data_heatmap[i][hue_order1[0]])
    data_heatmap_normalized.loc[0] = temp
    data_heatmap_normalized.index = [hue_order1[1]]
    data_heatmap_normalized.columns = ['%.3f' % elem for elem in x[:up]]
    df_int_normalized.loc[len(df_int_normalized.index)] = temp

    plt.subplots(figsize=(12, len(hue_order)))
    ax1 = sns.heatmap(data_heatmap_normalized, cbar=0, linewidths=2, vmax=data_heatmap_normalized.values.max(),
                      vmin=data_heatmap_normalized.values.min(), square=True, cmap='coolwarm')
    plt.savefig('%s/%s_heatmap_DNAFISH_int_n%s_cell%s_%s.pdf' % (output_dir, figure_name, n_dilation, len(data), hue_order[k]))
    plt.show()

# int_seg
for k in range(len(hue_order)):
    data = df[df['sample'] == hue_order[k]].copy().reset_index(drop=True)
    data_r = df_r[df_r['sample'] == hue_order[k]].copy().reset_index(drop=True)
    data_heatmap = pd.DataFrame(columns=x[:up])

    for s in hue_order1:
        data_sample = data_r[data_r['sample_category'] == s].copy().reset_index(drop=True)
        data_radial = []
        for i in range(len(x[:up])):
            if x[i] == 0.025:
                filter = (data_sample['relative_r'] >= 0) & (data_sample['relative_r'] <= 0.05)
            else:
                filter = (data_sample['relative_r'] > (x[i]-0.025)) & (data_sample['relative_r'] <= (x[i]+0.025))
            if s == hue_order1[0]:
                feature = 'seg'
            else:
                feature = 'intensity'
            data_radial.append(sum(data_sample[filter][feature])/sum(data_sample[feature]))
        data_heatmap.loc[len(data_heatmap.index)] = data_radial
    data_heatmap.index = hue_order1
    data_heatmap.columns = ['%.3f' % elem for elem in x[:up]]

    plt.subplots(figsize=(12, len(hue_order1)))
    ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
    plt.savefig('%s/%s_heatmap_DNAFISH_int-seg_before_normalization_n%s_cell%s_%s.pdf' % (output_dir, figure_name, n_dilation, len(data), hue_order[k]))
    plt.show()

    data_heatmap_normalized = pd.DataFrame(columns=x[:up])
    temp = []
    for i in data_heatmap.columns:
        temp.append(data_heatmap[i][hue_order1[1]] / data_heatmap[i][hue_order1[0]])
    data_heatmap_normalized.loc[0] = temp
    data_heatmap_normalized.index = [hue_order1[1]]
    data_heatmap_normalized.columns = ['%.3f' % elem for elem in x[:up]]
    df_int_seg_normalized.loc[len(df_int_seg_normalized.index)] = temp

    plt.subplots(figsize=(12, len(hue_order)))
    ax1 = sns.heatmap(data_heatmap_normalized, cbar=0, linewidths=2, vmax=data_heatmap_normalized.values.max(),
                      vmin=data_heatmap_normalized.values.min(), square=True, cmap='coolwarm')
    plt.savefig('%s/%s_heatmap_DNAFISH_int-seg_n%s_cell%s_%s.pdf' % (output_dir, figure_name, n_dilation, len(data), hue_order[k]))
    plt.show()

df_seg_normalized.index = hue_order
df_int_normalized.index = hue_order
df_int_seg_normalized.index = hue_order
df_seg_normalized.columns = ['%.3f' % elem for elem in x[:up]]
df_int_normalized.columns = ['%.3f' % elem for elem in x[:up]]
df_int_seg_normalized.columns = ['%.3f' % elem for elem in x[:up]]

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=df_seg_normalized.values.max(),
                  vmin=df_seg_normalized.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/%s_heatmap_DNAFISH_seg_n%s_cell%s.pdf' % (output_dir, figure_name, n_dilation, len(data)))
plt.show()

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(df_int_normalized, cbar=0, linewidths=2, vmax=df_int_normalized.values.max(),
                  vmin=df_int_normalized.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/%s_heatmap_DNAFISH_int_n%s_cell%s.pdf' % (output_dir, figure_name, n_dilation, len(data)))
plt.show()

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(df_int_seg_normalized, cbar=0, linewidths=2, vmax=df_int_seg_normalized.values.max(),
                  vmin=df_int_seg_normalized.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/%s_heatmap_DNAFISH_int-seg_n%s_cell%s.pdf' % (output_dir, figure_name, n_dilation, len(data)))
plt.show()

plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['sample'] == hue_order[k]].copy().reset_index(drop=True)
    plt.plot(x, df_seg_normalized.iloc[k].tolist(), color=line_colors[k], label='%s, n=%s' % (hue_order[k], len(data)))
plt.axhline(y=1, linestyle='--', color='r')
plt.xlabel('relative_r')
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s_rc-new_DNAFISH_seg_n%s.pdf' % (output_dir, figure_name, n_dilation))
plt.show()

plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['sample'] == hue_order[k]].copy().reset_index(drop=True)
    plt.plot(x, df_int_normalized.iloc[k].tolist(), color=line_colors[k], label='%s, n=%s' % (hue_order[k], len(data)))
plt.axhline(y=1, linestyle='--', color='r')
plt.xlabel('relative_r')
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s_rc-new_DNAFISH_int_n%s.pdf' % (output_dir, figure_name, n_dilation))
plt.show()

plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['sample'] == hue_order[k]].copy().reset_index(drop=True)
    plt.plot(x, df_int_seg_normalized.iloc[k].tolist(), color=line_colors[k], label='%s, n=%s' % (hue_order[k], len(data)))
plt.axhline(y=1, linestyle='--', color='r')
plt.xlabel('relative_r')
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s_rc-new_DNAFISH_int-seg_n%s.pdf' % (output_dir, figure_name, n_dilation))
plt.show()