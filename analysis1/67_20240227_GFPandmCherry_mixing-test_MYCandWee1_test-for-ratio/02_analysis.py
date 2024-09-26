import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240227_analysis_GFPandmCherry_mixing-ratio_test/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sprocessed/" % master_folder

folder = '384well'
sample = 'XY01'
fov = 15

data = pd.read_csv('%s/%s/txt/%s_%s_manual.txt' % (data_dir, folder, sample, fov), na_values=['.'], sep='\t')
data['log10_GFP'] = np.log10(data['GFP'])
data['log10_mCherry'] = np.log10(data['mCherry'])

print(len(data))
print(len(data[data['log10_mCherry']>3.25]))

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=data, x='log10_mCherry', y='log10_GFP', alpha=0.8, s=10)
# plt.xlim([0, 3])
# plt.ylim([0, 5000])
# plt.savefig('%s/%s_%s_MK1775_update1_hc.pdf' % (output_dir, folder, feature))
plt.show()
