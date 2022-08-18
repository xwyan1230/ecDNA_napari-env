import matplotlib.pyplot as plt
import pandas as pd
import shared.math as mat
import numpy as np
import seaborn as sns

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220809_Aarohi_chromosomeRecombination_metaphase" \
                "/08052022_ChrRecomb_KCl_timept/"
save_folder = master_folder
hist_colors = [(0.70, 0.70, 0.70), (0.65, 0.90, 0.90), (0.95, 0.50, 0.50)]
line_colors = [(0.50, 0.50, 0.50), (0.30, 0.70, 0.70), (0.85, 0.35, 0.25)]
sns.set_palette(sns.color_palette(hist_colors))

sample = 'ATCC'
sample1 = 'Jun16'
sample2 = 'XY'

data = pd.read_csv(("%s%s.txt" % (master_folder, sample)), na_values=['.'], sep='\t')
data1 = pd.read_csv(("%s%s.txt" % (master_folder, sample1)), na_values=['.'], sep='\t')
data2 = pd.read_csv(("%s%s.txt" % (master_folder, sample2)), na_values=['.'], sep='\t')
data['sample'] = [sample] * len(data)
data1['sample'] = [sample1] * len(data1)
data2['sample'] = [sample2] * len(data2)

ax1 = plt.scatter(data['DM_copy'], data['DM_total_int'], s=5, label=sample, color=hist_colors[0])
_, int_fit_r2_sample, int_fit_a_sample = mat.fitting_linear_b0(data['DM_copy'].tolist(), data['DM_total_int'].tolist())
plt.scatter(data1['DM_copy'], data1['DM_total_int'], s=5, label=sample1, color=hist_colors[1])
_, int_fit_r2_sample1, int_fit_a_sample1 = mat.fitting_linear_b0(data1['DM_copy'].tolist(), data1['DM_total_int'].tolist())
plt.scatter(data2['DM_copy'], data2['DM_total_int'], s=5, label=sample2, color=hist_colors[2])
_, int_fit_r2_sample2, int_fit_a_sample2 = mat.fitting_linear_b0(data2['DM_copy'].tolist(), data2['DM_total_int'].tolist())

x = np.arange(0, 200, 10)
y_sample = int_fit_a_sample * x
y_sample1 = int_fit_a_sample1 * x
y_sample2 = int_fit_a_sample2 * x
plt.plot(x, y_sample, linewidth=2, linestyle='--', color=line_colors[0])
plt.plot(x, y_sample1, linewidth=2, linestyle='--', color=line_colors[1])
plt.plot(x, y_sample2, linewidth=2, linestyle='--', color=line_colors[2])

plt.xlabel('DM copy number')
plt.ylabel('DM intensity')
plt.legend()
plt.savefig('%s/linear.pdf' % save_folder)
plt.close()

data['HSR_corrected_copy_number'] = data['HSR_total_int']/int_fit_a_sample[0]
data1['HSR_corrected_copy_number'] = data1['HSR_total_int']/int_fit_a_sample1[0]
data2['HSR_corrected_copy_number'] = data2['HSR_total_int']/int_fit_a_sample2[0]

plt.scatter(data['DM_copy'], data['HSR_corrected_copy_number'], s=5, label=sample, color=hist_colors[0])
plt.scatter(data1['DM_copy'], data1['HSR_corrected_copy_number'], s=5, label=sample1, color=hist_colors[1])
plt.scatter(data2['DM_copy'], data2['HSR_corrected_copy_number'], s=5, label=sample2, color=hist_colors[2])
plt.axhline(10, linestyle='--', color='r')
plt.axvline(5, linestyle='--', color='r')
plt.xlabel('DM copy number')
plt.ylabel('HSR copy number')
plt.legend()
plt.savefig('%s/scatter.pdf' % save_folder)
plt.close()

data_DM = data[(data['DM_copy'] > 5) & (data['HSR_corrected_copy_number'] < 10)]
data_DM1 = data1[(data1['DM_copy'] > 5) & (data1['HSR_corrected_copy_number'] < 10)]
data_DM2 = data2[(data2['DM_copy'] > 5) & (data2['HSR_corrected_copy_number'] < 10)]

data_sum = pd.concat([data_DM, data_DM1, data_DM2], axis=0, ignore_index=True)

histrange = [0, 300]
plt.hist(data_DM1['DM_copy'], bins=30, range=histrange, label=sample1, color=hist_colors[1], alpha=0.7)
plt.hist(data_DM2['DM_copy'], bins=30, range=histrange, label=sample2, color=hist_colors[2], alpha=0.7)
plt.hist(data_DM['DM_copy'], bins=30, range=histrange, label=sample, color=hist_colors[0], alpha=0.7)
plt.xlabel("copy number")
plt.legend()
plt.savefig('%s/DM_copy.pdf' % save_folder)
plt.close()

sns.set_palette(sns.color_palette(hist_colors))
sns.displot(data_sum, x="DM_copy", hue="sample", multiple="dodge")
plt.savefig('%s/DM_copy_dodge.pdf' % save_folder)
plt.close()

sns.set_palette(sns.color_palette(line_colors))
sns.displot(data_sum, x="DM_copy", hue="sample", kind='kde', fill='True')
plt.savefig('%s/DM_copy_kde.pdf' % save_folder)
plt.close()

print("DONE!")
