import matplotlib.pyplot as plt
import pandas as pd
import shared.math as mat
import shared.dataframe as dat
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
data['total_copy'] = data['DM_copy'] + data['HSR_corrected_copy_number']
data1['total_copy'] = data1['DM_copy'] + data1['HSR_corrected_copy_number']
data2['total_copy'] = data2['DM_copy'] + data2['HSR_corrected_copy_number']

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

data_sum = pd.concat([data, data1, data2], axis=0, ignore_index=True)

sns.set_palette(sns.color_palette(hist_colors))
sns.displot(data_sum, x="total_copy", hue="sample", multiple="dodge")
plt.savefig('%s/total_copy_dodge.pdf' % save_folder)
plt.close()

sns.set_palette(sns.color_palette(line_colors))
sns.displot(data_sum, x="total_copy", hue="sample", kind='kde', fill='True')
plt.savefig('%s/total_copy_kde.pdf' % save_folder)
plt.close()

print(sample)
mean_copy = np.sum(data['total_copy'].tolist())/len(data)
print(mean_copy)
print(sample1)
mean_copy1 = np.sum(data1['total_copy'].tolist())/len(data1)
print(mean_copy1)
print(sample2)
mean_copy2 = np.sum(data2['total_copy'].tolist())/len(data2)
print(mean_copy2)

sample_qPCR_copy = 376 * 2
sample1_qPCR_copy = 262.9 * 2

_, r2_sample, a_sample = mat.fitting_linear_b0([mean_copy, mean_copy1], [sample_qPCR_copy, sample1_qPCR_copy])
x = np.arange(0, 150, 10)
y = a_sample * x
print(a_sample)

plt.scatter([mean_copy, mean_copy1], [sample_qPCR_copy, sample1_qPCR_copy])
plt.plot(x, y, linewidth=2, linestyle='--', color='r')
plt.xlim([0, 150])
plt.ylim([0, 800])
plt.xlabel('copy_number_from_metaphase_spread')
plt.ylabel('copy_number_from_qPCR')
plt.savefig('%s/copy_number_comparison.pdf' % save_folder)
plt.close()

data_DM = data[(data['DM_copy'] > 5) & (data['HSR_corrected_copy_number'] < 10)]
data_DM1 = data1[(data1['DM_copy'] > 5) & (data1['HSR_corrected_copy_number'] < 10)]
data_DM2 = data2[(data2['DM_copy'] > 5) & (data2['HSR_corrected_copy_number'] < 10)]

data_DM_sum = pd.concat([data_DM, data_DM1, data_DM2], axis=0, ignore_index=True)

histrange = [0, 300]
plt.hist(data_DM1['DM_copy'], bins=30, range=histrange, label=sample1, color=hist_colors[1], alpha=0.7)
plt.hist(data_DM2['DM_copy'], bins=30, range=histrange, label=sample2, color=hist_colors[2], alpha=0.7)
plt.hist(data_DM['DM_copy'], bins=30, range=histrange, label=sample, color=hist_colors[0], alpha=0.7)
plt.xlabel("copy number")
plt.legend()
plt.savefig('%s/DM_copy.pdf' % save_folder)
plt.close()

sns.set_palette(sns.color_palette(hist_colors))
sns.displot(data_DM_sum, x="DM_copy", hue="sample", multiple="dodge")
plt.savefig('%s/DM_copy_dodge.pdf' % save_folder)
plt.close()

sns.set_palette(sns.color_palette(line_colors))
sns.displot(data_DM_sum, x="DM_copy", hue="sample", kind='kde', fill='True')
plt.savefig('%s/DM_copy_kde.pdf' % save_folder)
plt.close()

features = ['DM_ind_area', 'DM_ind_mean_int', 'DM_ind_total_int', 'HSR_ind_area', 'HSR_ind_mean_int', 'HSR_ind_total_int']
dataframes = [data, data1, data2]
for f in features:
    for df in dataframes:
        df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

data_area_ind_DM = dat.list_addup_from_df(data, 'DM_ind_area')
data_mean_int_ind_DM = dat.list_addup_from_df(data, 'DM_ind_mean_int')
data_total_int_ind_DM = dat.list_addup_from_df(data, 'DM_ind_total_int')
data_area_ind_HSR = dat.list_addup_from_df(data, 'HSR_ind_area')
data_mean_int_ind_HSR = dat.list_addup_from_df(data, 'HSR_ind_mean_int')
data_total_int_ind_HSR = dat.list_addup_from_df(data, 'HSR_ind_total_int')
data1_area_ind_DM = dat.list_addup_from_df(data1, 'DM_ind_area')
data1_mean_int_ind_DM = dat.list_addup_from_df(data1, 'DM_ind_mean_int')
data1_total_int_ind_DM = dat.list_addup_from_df(data1, 'DM_ind_total_int')
data1_area_ind_HSR = dat.list_addup_from_df(data1, 'HSR_ind_area')
data1_mean_int_ind_HSR = dat.list_addup_from_df(data1, 'HSR_ind_mean_int')
data1_total_int_ind_HSR = dat.list_addup_from_df(data1, 'HSR_ind_total_int')
data2_area_ind_DM = dat.list_addup_from_df(data2, 'DM_ind_area')
data2_mean_int_ind_DM = dat.list_addup_from_df(data2, 'DM_ind_mean_int')
data2_total_int_ind_DM = dat.list_addup_from_df(data2, 'DM_ind_total_int')
data2_area_ind_HSR = dat.list_addup_from_df(data2, 'HSR_ind_area')
data2_mean_int_ind_HSR = dat.list_addup_from_df(data2, 'HSR_ind_mean_int')
data2_total_int_ind_HSR = dat.list_addup_from_df(data2, 'HSR_ind_total_int')

data_m_DM_sum = pd.DataFrame({'sample': [sample] * len(data_area_ind_DM) + [sample1] * len(data1_area_ind_DM) +
                                        [sample2] * len(data2_area_ind_DM),
                              'DM_ind_area': data_area_ind_DM + data1_area_ind_DM + data2_area_ind_DM,
                              'DM_ind_mean_int': data_mean_int_ind_DM + data1_mean_int_ind_DM + data2_mean_int_ind_DM,
                              'DM_ind_total_int': data_total_int_ind_DM + data1_total_int_ind_DM +
                                                  data2_total_int_ind_DM})
data_m_HSR_sum = pd.DataFrame({'sample': [sample] * len(data_area_ind_HSR) + [sample1] * len(data1_area_ind_HSR) +
                                         [sample2] * len(data2_area_ind_HSR),
                               'HSR_ind_area': data_area_ind_HSR + data1_area_ind_HSR + data2_area_ind_HSR,
                               'HSR_ind_mean_int': data_mean_int_ind_HSR + data1_mean_int_ind_HSR +
                                                   data2_mean_int_ind_HSR,
                               'HSR_ind_total_int': data_total_int_ind_HSR + data1_total_int_ind_HSR +
                                                   data2_total_int_ind_HSR})

DM_plot_features = ['DM_ind_area', 'DM_ind_mean_int', 'DM_ind_total_int']
HSR_plot_features = ['HSR_ind_area', 'HSR_ind_mean_int', 'HSR_ind_total_int']

for f in DM_plot_features:
    sns.set_palette(sns.color_palette(hist_colors))
    sns.displot(data_m_DM_sum, x=f, hue="sample", multiple="dodge")
    plt.savefig('%s/%s_dodge.pdf' % (save_folder, f))
    plt.close()

    sns.set_palette(sns.color_palette(line_colors))
    sns.displot(data_m_DM_sum, x=f, hue="sample", kind='kde', fill='True')
    plt.savefig('%s/%s_kde.pdf' % (save_folder, f))
    plt.close()

for f in HSR_plot_features:
    sns.set_palette(sns.color_palette(hist_colors))
    sns.displot(data_m_HSR_sum, x=f, hue="sample", multiple="dodge")
    plt.savefig('%s/%s_dodge.pdf' % (save_folder, f))
    plt.close()

    sns.set_palette(sns.color_palette(line_colors))
    sns.displot(data_m_HSR_sum, x=f, hue="sample", kind='kde', fill='True')
    plt.savefig('%s/%s_kde.pdf' % (save_folder, f))
    plt.close()

f = 'DM_ind_area'
sns.set_palette(sns.color_palette(hist_colors))
sns.displot(data_m_DM_sum, x=f, hue="sample", multiple="dodge")
plt.xlim([0, 200])
plt.savefig('%s/%s_dodge_xlim.pdf' % (save_folder, f))
plt.close()

sns.set_palette(sns.color_palette(line_colors))
sns.displot(data_m_DM_sum, x=f, hue="sample", kind='kde', fill='True')
plt.xlim([0, 200])
plt.savefig('%s/%s_kde_xlim.pdf' % (save_folder, f))
plt.close()



print("DONE!")
