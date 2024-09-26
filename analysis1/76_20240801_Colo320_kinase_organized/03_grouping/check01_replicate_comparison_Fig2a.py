# Fig 2a

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
# output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/figures/"
batch = 'point5uM_48hr'
feature = 'log2fc_hoechst'

# PARAMETERS
survival_features = ['log2fc_hoechst']
cc_features = []
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']

data = pd.read_csv('%s/01_summary/%s_grouping.txt' % (data_dir, batch), na_values=['.'], sep='\t')
if feature in survival_features:
    data_filter = 'survival_filter'
elif feature in cc_features:
    data_filter = 'cc_filter'
else:
    data_filter = 0

if data_filter != 0:
    data_flt = data[data[data_filter] == 1].copy().reset_index(drop=True)
else:
    data_flt = data.copy()

print(len(data_flt))

fig, ax = plt.subplots(figsize=(9, 7))
# fig.subplots_adjust(right=0.8)
sns.scatterplot(data=data_flt, x='%s_rep1' % feature, y='%s_rep2' % feature, s=20, alpha=1, color='#cccccc')
sns.scatterplot(data=data_flt[data_flt['treatment']=='DMSO'], x='%s_rep1' % feature, y='%s_rep2' % feature, s=20, alpha=1, color=rainboo_colors[6])
# plt.xlim([7, 16.5])
# plt.ylim([9, 16.5])
plt.show()

print("DONE!")