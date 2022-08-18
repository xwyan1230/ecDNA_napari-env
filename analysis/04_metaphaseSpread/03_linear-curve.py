import matplotlib.pyplot as plt
import pandas as pd
import shared.math as mat
import numpy as np
import seaborn as sns
import shared.dataframe as dat

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220809_Aarohi_chromosomeRecombination_metaphase" \
                "/08052022_ChrRecomb_KCl_timept/"
save_folder = master_folder

sample = 'ATCC'

data = pd.read_csv(("%s%s.txt" % (master_folder, sample)), na_values=['.'], sep='\t')

ax1 = plt.scatter(data['DM_copy'], data['DM_total_int'], s=5, label=sample)
_, int_fit_r2_sample, int_fit_a_sample = mat.fitting_linear_b0(data['DM_copy'].tolist(), data['DM_total_int'].tolist())

print(int_fit_r2_sample)
print(int_fit_a_sample[0])

x = np.arange(0, 200, 10)
y_sample = int_fit_a_sample * x
plt.plot(x, y_sample, linewidth=2, linestyle='--')

plt.xlabel('DM copy number')
plt.ylabel('DM intensity')
plt.legend()
plt.savefig('%s/%s_linear.pdf' % (save_folder, sample))
plt.close()

data['HSR_corrected_copy_number'] = data['HSR_total_int']/int_fit_a_sample[0]

plt.scatter(data['DM_copy'], data['HSR_corrected_copy_number'], s=5, label=sample)
plt.axhline(10, linestyle='--', color='r')
plt.axvline(5, linestyle='--', color='r')
plt.xlabel('DM copy number')
plt.ylabel('HSR copy number')
plt.legend()
plt.savefig('%s/%s_scatter.pdf' % (save_folder, sample))
plt.close()

data_DM = data[(data['DM_copy'] > 5) & (data['HSR_corrected_copy_number'] < 10)]
n_DM = len(data_DM)
n_HSR = len(data[(data['DM_copy'] < 5) & (data['HSR_corrected_copy_number'] > 10)])
n_hybrid = len(data[(data['DM_copy'] >= 5) & (data['HSR_corrected_copy_number'] >= 10)])
n_zero = len(data[(data['DM_copy'] <= 5) & (data['HSR_corrected_copy_number'] <= 10)])
n_total = len(data)

print(n_DM/n_total)
print(n_HSR/n_total)
print(n_hybrid/n_total)
print(n_zero/n_total)

plt.hist(data_DM['DM_copy'], bins=25)
plt.xlabel("copy number")
plt.legend()
plt.savefig('%s/%s_DM_copy.pdf' % (save_folder, sample))
plt.close()

sns.displot(data_DM, x="DM_copy", kde=True)
plt.savefig('%s/DM_copy_kde.pdf' % save_folder)
plt.close()

data_area_ind_DM = dat.list_addup_from_df(data, 'DM_ind_area')
data_mean_int_ind_DM = dat.list_addup_from_df(data, 'DM_ind_mean_int')
data_total_int_ind_DM = dat.list_addup_from_df(data, 'DM_ind_total_int')
data_area_ind_HSR = dat.list_addup_from_df(data, 'HSR_ind_area')
data_mean_int_ind_HSR = dat.list_addup_from_df(data, 'HSR_ind_mean_int')
data_total_int_ind_HSR = dat.list_addup_from_df(data, 'HSR_ind_total_int')

print("DONE!")
