import pandas as pd
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220325_Natasha_THZ1/"
sample = ['72hr_100nMTHZ', 'DMSO']
hist_colors = [(0.90, 0.90, 0.90), (0.95, 0.50, 0.50)]
line_colors = [(0.50, 0.50, 0.50), (0.85, 0.35, 0.25)]
feature = ['radial_curve_DNAFISH', 'radial_curve_nuclear', 'angle_curve_DNAFISH', 'angle_curve_nuclear',
           'area_individual_ecDNA', 'mean_int_individual_ecDNA',
           'percentage_area_curve_ecDNA', 'percentage_int_curve_ecDNA']

# LOAD FILE
data1 = pd.read_csv('%s%s/%s.txt' % (master_folder, sample[0], sample[0]), na_values=['.'], sep='\t')
data2 = pd.read_csv('%s%s/%s.txt' % (master_folder, sample[1], sample[1]), na_values=['.'], sep='\t')
for f in feature:
    data1[f] = [dat.str_to_float(data1[f][i]) for i in range(len(data1))]
    data2[f] = [dat.str_to_float(data2[f][i]) for i in range(len(data2))]

# calculate
data1['n_ecDNA'] = [len(data1['percentage_area_curve_ecDNA'][i])-1 for i in range(len(data1))]
data2['n_ecDNA'] = [len(data2['percentage_area_curve_ecDNA'][i])-1 for i in range(len(data2))]
data1_area_ind_ecDNA = dat.list_addup_from_df(data1, 'area_individual_ecDNA')
data2_area_ind_ecDNA = dat.list_addup_from_df(data2, 'area_individual_ecDNA')
data1['cum_area_ind_ecDNA'] = [dat.list_sum(data1['area_individual_ecDNA'][i]) for i in range(len(data1))]
data2['cum_area_ind_ecDNA'] = [dat.list_sum(data2['area_individual_ecDNA'][i]) for i in range(len(data2))]
max_n_ecDNA1 = max(data1['n_ecDNA'])
max_n_ecDNA2 = max(data2['n_ecDNA'])
data1['cum_area_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data1['cum_area_ind_ecDNA'])
data2['cum_area_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data2['cum_area_ind_ecDNA'])

# hist0
"""f = 'n_ecDNA'

plt.subplots(figsize=(6, 4))
weights1 = np.ones_like(data1[f]) / len(data1)
weights2 = np.ones_like(data2[f]) / len(data2)
plt.hist([data1[f], data2[f]], weights=[weights1, weights2], color=colors,
         edgecolor=(0.2, 0.2, 0.2), label=sample, bins=20)
plt.xlabel(f)
plt.ylabel('Probability')
plt.legend()
plt.show()"""

# hist1

"""d1 = data1_area_ind_ecDNA
d2 = data2_area_ind_ecDNA
xlabel = 'area_individual_ecDNA'

plt.subplots(figsize=(6, 4))
weights1 = np.ones_like(d1) / len(d1)
weights2 = np.ones_like(d2) / len(d2)
plt.hist([d1, d2], weights=[weights1, weights2], color=colors,
         edgecolor=(0.2, 0.2, 0.2), label=sample, bins=20)
plt.xlabel(xlabel)
plt.ylabel('Probability')
plt.ylim([0, 0.1])
plt.legend()
plt.show()"""

# curve0
"""f = 'radial_curve_DNAFISH'
fc = 'radial_curve_nuclear'
# fc = 0
x = np.arange(0, 0.98, 0.01)
x_label = 'r'

number_nuclear1 = len(data1)
number_nuclear2 = len(data2)

mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data1[f].tolist())
mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(data2[f].tolist())
if fc != 0:
    mean_curve1c, ci_lower1c, ci_higher1c = dat.mean_list(data1[fc].tolist())
    mean_curve2c, ci_lower2c, ci_higher2c = dat.mean_list(data2[fc].tolist())

plt.subplots(figsize=(6, 4))
for i in range(len(data1)):
    plt.plot(x, data1[f][i], alpha=0.05, color=[colors[0][j]+0.1 for j in range(len(colors[0]))])
for i in range(len(data2)):
    plt.plot(x, data2[f][i], alpha=0.05, color=[colors[1][j]+0.1 for j in range(len(colors[1]))])
plt.plot(x, mean_curve1, color=colors[0], label='%s, n=%s' % (sample[0], number_nuclear1))
plt.plot(x, mean_curve2, color=colors[1], label='%s, n=%s' % (sample[1], number_nuclear2))
plt.plot(x, ci_lower1, color=colors[0], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher1, color=colors[0], linestyle='--', linewidth=0.5)
plt.plot(x, ci_lower2, color=colors[1], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher2, color=colors[1], linestyle='--', linewidth=0.5)
if fc != 0:
    plt.plot(x, mean_curve1c, color=[colors[0][j]+0.15 for j in range(len(colors[0]))])
    plt.plot(x, mean_curve2c, color=[colors[1][j]+0.15 for j in range(len(colors[1]))])
plt.xlabel(x_label)
plt.ylabel(f)
plt.legend()
plt.savefig('%s/radial_comparison_%s_vs_%s.pdf' % (master_folder, sample[0], sample[1]))
plt.close()"""

# curve1
"""f = 'angle_curve_DNAFISH'
# fc = 'angle_curve_nuclear'
fc = 0
x = np.arange(0, 360, 1)
x_label = 'degree'

number_nuclear1 = len(data1)
number_nuclear2 = len(data2)

mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data1[f].tolist())
mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(data2[f].tolist())
if fc != 0:
    mean_curve1c, ci_lower1c, ci_higher1c = dat.mean_list(data1[fc].tolist())
    mean_curve2c, ci_lower2c, ci_higher2c = dat.mean_list(data2[fc].tolist())

plt.subplots(figsize=(6, 4))
for i in range(len(data1)):
    plt.plot(x, data1[f][i], alpha=0.05, color=[colors[0][j]+0.1 for j in range(len(colors[0]))])
for i in range(len(data2)):
    plt.plot(x, data2[f][i], alpha=0.05, color=[colors[1][j]+0.1 for j in range(len(colors[1]))])
plt.plot(x, mean_curve1, color=colors[0], label='%s, n=%s' % (sample[0], number_nuclear1))
plt.plot(x, mean_curve2, color=colors[1], label='%s, n=%s' % (sample[1], number_nuclear2))
plt.plot(x, ci_lower1, color=colors[0], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher1, color=colors[0], linestyle='--', linewidth=0.5)
plt.plot(x, ci_lower2, color=colors[1], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher2, color=colors[1], linestyle='--', linewidth=0.5)
if fc != 0:
    plt.plot(x, mean_curve1c, color=[colors[0][j]+0.15 for j in range(len(colors[0]))])
    plt.plot(x, mean_curve2c, color=[colors[1][j]+0.15 for j in range(len(colors[1]))])
plt.xlabel(x_label)
plt.ylabel(f)
plt.legend()
plt.savefig('%s/angle_comparison_%s_vs_%s.pdf' % (master_folder, sample[0], sample[1]))
plt.close()"""

# curve2
"""f = 'cum_area_ind_ecDNA_filled'
x_label = 'number of ecDNA hub'

# data1 = data1[data1['total_area_ecDNA'] > 100]
# data2 = data2[data2['total_area_ecDNA'] > 100]
# data1.reset_index(drop=True, inplace=True)
# data2.reset_index(drop=True, inplace=True)

number_nuclear1 = len(data1)
number_nuclear2 = len(data2)

mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(dat.list_fill_with_num(data1[f].tolist(), 1))
mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(dat.list_fill_with_num(data2[f].tolist(), 1))

plt.subplots(figsize=(6, 4))
for i in range(len(data1)):
    plt.plot(data1[f][i], alpha=0.05, color=[colors[0][j]+0.1 for j in range(len(colors[0]))])
for i in range(len(data2)):
    plt.plot(data2[f][i], alpha=0.05, color=[colors[1][j]+0.1 for j in range(len(colors[1]))])
plt.plot(mean_curve1, color=colors[0], label='%s, n=%s' % (sample[0], number_nuclear1))
plt.plot(mean_curve2, color=colors[1], label='%s, n=%s' % (sample[1], number_nuclear2))
plt.plot(ci_lower1, color=colors[0], linestyle='--', linewidth=0.5)
plt.plot(ci_higher1, color=colors[0], linestyle='--', linewidth=0.5)
plt.plot(ci_lower2, color=colors[1], linestyle='--', linewidth=0.5)
plt.plot(ci_higher2, color=colors[1], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylabel(f)
plt.legend()
plt.show()
# plt.savefig('%s/cum_area_comparison_%s_vs_%s.pdf' % (master_folder, sample[0], sample[1]))
# plt.close()"""
