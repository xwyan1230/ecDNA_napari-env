import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
data_dir1 = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/processed/"
# output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/figures/"
batch = '5uM_48hr'
hit_list = ['Ceritinib', 'Flavopiridol', 'Bosutinib', 'Ponatinib', 'Abemaciclib', 'Crenolanib', 'Torin 2', 'NVP-2',
            'PND-1186', 'GSK1838705A', 'ASK1-IN-1', 'Defactinib', 'GNF-7', 'AZ191', 'SGI-1776']

"""'Rapamycin', 'AZD-7762',
            'TL02-59', 'SCH900776', 'WNK463', 'NVP-AEW541', 'Ripasudil', 'Mps1-IN-3', 'URMC-099', 'CRT 0066101',
            'Mps1-IN-1', 'Mps1-IN-2', 'Bemcentinib', 'Upadacitinib', 'HPK1-IN-3', 'A-674563', 'XMU-MP-1', 'KU-60019'"""

# hit_list = ['Flavopiridol', 'NVP-2', 'BI 2536', 'GSK461364']
# hit_list = ['AZD-7762']
# hit_list = ['Mps1-IN-3', 'Mps1-IN-1', 'Mps1-IN-2']

# PARAMETERS
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']

data = pd.read_csv('%s/01_summary/%s_normalize.txt' % (data_dir, batch), na_values=['.'], sep='\t')
data = data[~data['treatment'].isin(['ctrl'])].copy().reset_index(drop=True)

plt.subplots(figsize=(9, 7))
feature = 'hoechst_neg_G1_normalized'
feature1 = 'hoechst_neg_G2M_normalized'
sns.scatterplot(data=data, x=feature, y=feature1, alpha=1, s=20, color='#cccccc')
sns.scatterplot(data=data[data['treatment'] == 'DMSO'], x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[6])
# data_flt = data[data['%s_norm' % feature1] < 1.35].copy().reset_index(drop=True)
data_flt = data[data['treatment'].isin(hit_list)].copy().reset_index(drop=True)
sns.scatterplot(data=data_flt, x=feature, y=feature1, alpha=1, s=20, color=rainboo_colors[2])
data_flt = data[data[feature] > 22500].copy().reset_index(drop=True)
for i in range(len(data_flt)):
    plt.text(x=data_flt[feature][i] + 0.1, y=data_flt[feature1][i], s=data_flt['treatment'][i],
             size=6, color=rainboo_colors[2])
# plt.ylim(ylim_val)
# plt.xlim(xlim_val)
plt.show()