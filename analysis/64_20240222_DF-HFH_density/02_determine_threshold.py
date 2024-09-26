import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240222_analysis_DF-HFH_density/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sprocessed/" % master_folder

folder = '48hr_384well'

sample_name = 'HFH-Ctrl_update'
data = pd.read_csv('%s/%s/txt/%s.txt' % (data_dir, folder, sample_name), na_values=['.'], sep='\t')
data['log10_emiRFP670'] = np.log10(data['emiRFP670'])

hc = 0
print(np.mean(data['hoechst']))
fig, ax = plt.subplots(figsize=(9, 6))
fig.subplots_adjust(right=0.8)
sns.histplot(data=data, x='hoechst')
plt.axvline(hc, 0, 1000, c='r')
# plt.legend(loc=(1.04, 0))
if not os.path.exists("%s/%s/hoechst_hist/" % (output_dir, folder)):
    os.makedirs("%s/%s/hoechst_hist/" % (output_dir, folder))
# plt.savefig('%s/%s/hoechst_hist/%s.pdf' % (output_dir, folder, sample_name))
plt.show()

lc = 3.1
uc = 3.1
print(len(data[data['hoechst'] > hc]))
print(len(data[(data['hoechst'] > hc) & (data['log10_emiRFP670'] < lc)]))
print(len(data[(data['hoechst'] > hc) & (data['log10_emiRFP670'] > uc)]))
print(len(data[(data['hoechst'] > hc) & (data['log10_emiRFP670'] < lc)]) / len(data[data['hoechst'] > hc]))
print(len(data[(data['hoechst'] > hc) & (data['log10_emiRFP670'] > uc)]) / len(data[data['hoechst'] > hc]))

fig, ax = plt.subplots(figsize=(9, 6))
sns.histplot(data=data[data['hoechst'] > hc], x='log10_emiRFP670', bins=20)  # annot=True
plt.xlim([2, 6])
plt.axvline(lc, 0, 1000, c='r')
plt.axvline(uc, 0, 1000, c='r')
# plt.axvline(16000, 0, 1000, c='r')
# plt.legend(loc=(1.04, 0))
if not os.path.exists("%s/%s/emiRFP670_hist/" % (output_dir, folder)):
    os.makedirs("%s/%s/emiRFP670_hist/" % (output_dir, folder))
# plt.savefig('%s/%s/emiRFP670_hist/%s.pdf' % (output_dir, folder, sample_name))
plt.show()

