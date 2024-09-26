# Fig 1e

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
data_dir1 = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/processed/"
# output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/figures/"
batches = ['24hr_density', '48hr_density']
mode = 'mix'  # only accepts 'neg', 'pos', 'mix'

# PARAMETERS
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
fov_area = 2.336  # mm2
well_area = 11
neg = ['XY%s' % (201 + 2*x) for x in range(5)]
pos = ['XY%s' % (202 + 2*x) for x in range(5)]

df = pd.DataFrame(columns=['batch', 'identity', 'group', 'percentage'])
for batch in batches:
    if batch == '24hr_density':
        density = 12
    else:
        density = 7

    data = pd.read_csv('%s/01_summary/%s_normalize.txt' % (data_dir, batch), na_values=['.'], sep='\t')
    if mode == 'neg':
        data_cc = data[data['sample'].isin(neg)].copy().reset_index(drop=True)
    elif mode == 'pos':
        data_cc = data[data['sample'].isin(pos)].copy().reset_index(drop=True)
    else:
        data_cc = data[data['density'] == density].copy().reset_index(drop=True)

    for i in range(len(data_cc)):
        df.loc[len(df.index)] = [batch, 'G1', '%s_G1' % batch, data_cc['n_G1'][i]/data_cc['n_filtered'][i]]
        df.loc[len(df.index)] = [batch, 'G1S', '%s_G1S' % batch, data_cc['n_G1S'][i] / data_cc['n_filtered'][i]]
        df.loc[len(df.index)] = [batch, 'S', '%s_S' % batch, data_cc['n_S'][i] / data_cc['n_filtered'][i]]
        df.loc[len(df.index)] = [batch, 'G2M', '%s_G2M' % batch, data_cc['n_G2M'][i] / data_cc['n_filtered'][i]]
        df.loc[len(df.index)] = [batch, 'G2MG1', '%s_G2MG1' % batch, data_cc['n_G2MG1'][i] / data_cc['n_filtered'][i]]

plt.subplots(figsize=(9, 7))
sns.barplot(data=df, x='group', y='percentage', color='#e3e3e3', edgecolor='#4d4e4e')
sns.swarmplot(data=df, x='group', y='percentage', color='#bc4d4a', s=5)
plt.ylim([0, 0.5])
plt.show()