import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240209_analysis_Colo320_kinaseScreen_point5uM/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

folder = 'R2P1'

data = pd.read_csv('%s/%s/txt/XY07.txt' % (data_dir, folder), na_values=['.'], sep='\t')
data['log10_emiRFP670'] = np.log10(data['emiRFP670'])

hc = 11000
print(np.mean(data['hoechst']))
print(np.mean(data[data['hoechst']>hc]['hoechst']))
fig, ax = plt.subplots(figsize=(9, 6))
fig.subplots_adjust(right=0.8)
sns.histplot(data=data, x='hoechst')
plt.axvline(hc, 0, 1000, c='r')
# plt.legend(loc=(1.04, 0))
# plt.savefig('%s/%s/hoechst_hist/XY%s.pdf' % (output_dir, folder, i+1))
plt.show()

lc = 2.95
uc = 2.95
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
# plt.savefig('%s/%s/hoechst_hist/XY%s.pdf' % (output_dir, folder, i+1))
plt.show()

