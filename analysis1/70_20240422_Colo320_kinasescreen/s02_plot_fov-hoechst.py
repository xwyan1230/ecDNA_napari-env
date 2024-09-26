import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import seaborn as sns
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

if not os.path.exists("%sfigures/" % master_folder):
    os.makedirs("%sfigures/" % master_folder)

batch = '5uM_48hr'
plates = [1, 2]
rep = 1
cell = 'DF+HFH'

# target = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % data_dir, na_values=['.'])
# print(target.head())

data = pd.DataFrame()
for plate in plates:
    df = pd.read_csv('%s/%s/%s_%s.txt' % (data_dir, batch, batch, plate), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0)
data_screen = data[data['cell'] == cell].copy().reset_index(drop=True)
data_screen['log10_fov_hoechst'] = np.log10(data_screen['fov_hoechst'])
print(min(data_screen['log10_fov_hoechst']))
print(max(data_screen['log10_fov_hoechst']))

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_screen, x='n_neg', y='n_pos', c=data_screen['log10_fov_hoechst'], alpha=0.5, s=30, vmin=8.5, vmax=10.5)
plt.xlim([0, 2200])
plt.ylim([0, 2200])
plt.legend()
plt.savefig('%s/R%s_%s_n-neg_vs_n-pos_fov-hoechst.pdf' % (output_dir, rep, cell))
plt.show()

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_screen[data_screen['treatment'] == 'DMSO'], x='n_neg', y='n_pos', c=data_screen[data_screen['treatment'] == 'DMSO']['log10_fov_hoechst'], alpha=0.5, s=30, vmin=8.5, vmax=10.5)
plt.xlim([0, 2200])
plt.ylim([0, 2200])
plt.legend()
plt.savefig('%s/R%s_%s_n-neg_vs_n-pos_fov-hoechst_DMSO.pdf' % (output_dir, rep, cell))
plt.show()








