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

batch = 'point5uM_48hr'
exclude = ['Barasertib', 'Alisertib', 'AZD2811', 'PF-03814735', 'GSK269962A', 'Y-33075', 'BDP5290', 'Midostaurin', 'FF-10101']

data = pd.read_csv('%s/%s/%s_average.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data1 = pd.read_csv('%s/5uM_48hr/5uM_48hr_average.txt' % data_dir, na_values=['.'], sep='\t')
data['ref'] = data1['mean_hoechst']
data_flt = data[~data['treatment'].isin(exclude)].copy().reset_index(drop=True)
feature = 'mean_n'
mean_ctrl = np.mean(data_flt[data_flt['treatment']=='DMSO'][feature])
std_ctrl = np.std(data_flt[data_flt['treatment']=='DMSO'][feature])
data_flt = data_flt.sort_values(by='ref').reset_index(drop=True)

data_hp = pd.DataFrame()
data_hp[feature] = data_flt[feature]
data_hp.index = data_flt['treatment']

fig, ax = plt.subplots(figsize=(7, 60))
fig.subplots_adjust(left=0.3)
ax1 = sns.heatmap(data_hp, cbar=0, linewidths=2, vmax=mean_ctrl+10*std_ctrl, vmin=mean_ctrl-10*std_ctrl, square=True, cmap='coolwarm', annot=False, fmt='.2f')
plt.savefig('%s/%s/%s_heatmap_%s.pdf' % (output_dir, batch, batch, feature))
plt.show()