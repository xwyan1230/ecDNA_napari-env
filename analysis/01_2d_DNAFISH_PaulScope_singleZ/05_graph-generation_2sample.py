import pandas as pd
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt
import random
import seaborn as sns

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220325_Natasha_THZ1/"
sample = ['DMSO', '72hr_100nMTHZ']
hist_colors = [(0.90, 0.90, 0.90), (0.95, 0.50, 0.50)]
line_colors = [(0.50, 0.50, 0.50), (0.85, 0.35, 0.25)]
# '#FFA500', '#40E0D0', '#DA70D6', '#00CED1'
feature = ['radial_curve_DNAFISH', 'radial_curve_nuclear', 'angle_curve_DNAFISH', 'angle_curve_nuclear',
           'area_individual_ecDNA', 'mean_int_individual_ecDNA', 'total_int_individual_ecDNA',
           'percentage_area_curve_ecDNA', 'percentage_int_curve_ecDNA',
           'cum_area_ind_ecDNA', 'cum_area_ind_ecDNA_filled', 'cum_int_ind_ecDNA', 'cum_int_ind_ecDNA_filled']

# LOAD FILE
data1 = pd.read_csv('%s%s/%s.txt' % (master_folder, sample[0], sample[0]), na_values=['.'], sep='\t')
data2 = pd.read_csv('%s%s/%s.txt' % (master_folder, sample[1], sample[1]), na_values=['.'], sep='\t')
for f in feature:
    data1[f] = [dat.str_to_float(data1[f][i]) for i in range(len(data1))]
    data2[f] = [dat.str_to_float(data2[f][i]) for i in range(len(data2))]

# calculate
data1['area_individual_ecDNA_hub'] = dat.filter_small_from_lst_in_df(data1, 'area_individual_ecDNA', 100)
data2['area_individual_ecDNA_hub'] = dat.filter_small_from_lst_in_df(data2, 'area_individual_ecDNA', 100)
data1_area_ind_ecDNA = dat.list_addup_from_df(data1, 'area_individual_ecDNA')
data2_area_ind_ecDNA = dat.list_addup_from_df(data2, 'area_individual_ecDNA')

# hist0-batch
feature = ['n_ecDNA', 'max_area_ecDNA', 'area_nuclear', 'total_int_DNAFISH', 'total_int_nuclear', 'total_area_ecDNA',
           'total_int_ecDNA']

for i in feature:
    f = i
    plt.subplots(figsize=(6, 4))
    weights1 = np.ones_like(data1[f]) / len(data1)
    weights2 = np.ones_like(data2[f]) / len(data2)
    plt.hist([data1[f], data2[f]], weights=[weights1, weights2], color=hist_colors,
             edgecolor=(0.2, 0.2, 0.2), label=sample, bins=20)
    plt.xlabel(f)
    plt.ylabel('Probability')
    plt.legend()
    plt.savefig('%s/%s_%s_vs_%s.pdf' % (master_folder, f, sample[0], sample[1]))
    plt.close()

# hist1
d1 = [i for i in data1_area_ind_ecDNA if i > 100]
d2 = [i for i in data2_area_ind_ecDNA if i > 100]

xlabel = 'area_individual_ecDNA (>100)'

plt.subplots(figsize=(6, 4))
weights1 = np.ones_like(d1) / len(d1)
weights2 = np.ones_like(d2) / len(d2)
plt.hist([d1, d2], weights=[weights1, weights2], color=hist_colors,
         edgecolor=(0.2, 0.2, 0.2), label=sample, bins=20)
plt.xlabel(xlabel)
plt.ylabel('Probability')
plt.legend()
plt.savefig('%s/individual_ecDNA_%s_vs_%s.pdf' % (master_folder, sample[0], sample[1]))
plt.close()

# curve1
f = 'angle_curve_DNAFISH'
x = np.arange(0, 360, 1)
x_label = 'degree'

number_nuclear1 = len(data1)
number_nuclear2 = len(data2)

mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data1[f].tolist())
mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(data2[f].tolist())

plt.subplots(figsize=(6, 4))
for i in range(len(data1)):
    plt.plot(x, data1[f][i], alpha=0.05, color=[line_colors[0][j]+0.05 for j in range(len(line_colors[0]))])
for i in range(len(data2)):
    plt.plot(x, data2[f][i], alpha=0.05, color=[line_colors[1][j]+0.05 for j in range(len(line_colors[1]))])
plt.plot(x, mean_curve1, color=line_colors[0], label='%s, n=%s' % (sample[0], number_nuclear1))
plt.plot(x, mean_curve2, color=line_colors[1], label='%s, n=%s' % (sample[1], number_nuclear2))
plt.plot(x, ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x, ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylabel(f)
plt.legend()
plt.savefig('%s/%s_%s_vs_%s.pdf' % (master_folder, f, sample[0], sample[1]))
plt.close()

# curve2
feature = ['cum_int_ind_ecDNA_filled', 'cum_area_ind_ecDNA_filled']
for f in feature:
    x_label = 'number of ecDNA hub'

    number_nuclear1 = len(data1)
    number_nuclear2 = len(data2)

    mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(dat.list_fill_with_num(data1[f].tolist(), 1))
    mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(dat.list_fill_with_num(data2[f].tolist(), 1))

    plt.subplots(figsize=(6, 4))
    for i in range(len(data1)):
        plt.plot(data1[f][i], alpha=0.05, color=[line_colors[0][j]+0.05 for j in range(len(line_colors[0]))])
    for i in range(len(data2)):
        plt.plot(data2[f][i], alpha=0.05, color=[line_colors[1][j]+0.05 for j in range(len(line_colors[1]))])
    plt.plot(mean_curve1, color=line_colors[0], label='%s, n=%s' % (sample[0], number_nuclear1))
    plt.plot(mean_curve2, color=line_colors[1], label='%s, n=%s' % (sample[1], number_nuclear2))
    plt.plot(ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.xlabel(x_label)
    plt.ylabel(f)
    plt.legend()
    # plt.show()
    plt.savefig('%s/%s_%s_vs_%s.pdf' % (master_folder, f, sample[0], sample[1]))
    plt.close()
