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

batch = 'point5uM_48hr'
plates = [1, 2, 3]
rep = [1, 2]
cell = 'DF+HFH'

# target = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % data_dir, na_values=['.'])
# print(target.head())

data = pd.DataFrame()
for plate in plates:
    df = pd.read_csv('%s/%s/%s_%s.txt' % (data_dir, batch, batch, plate), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0)
data_screen = data[(data['cell'] == cell) & (data['rep'].isin(rep))].copy().reset_index(drop=True)
data_ctrl = data_screen[data_screen['treatment'] == 'DMSO'].copy().reset_index(drop=True)

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_screen, x='n_filtered', y='fov_hoechst', alpha=0.5, s=30)
sns.scatterplot(data=data_ctrl, x='n_filtered', y='fov_hoechst', color='r', alpha=0.5, s=30)
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_n-filtered_vs_fov-hoechst_individual.pdf' % (output_dir, batch, batch))
plt.show()








