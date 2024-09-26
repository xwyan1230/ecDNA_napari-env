# Fig 2g

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import utilities as uti
from scipy.stats import pearsonr

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
batches = ['point5uM_24hr', 'point5uM_48hr', '5uM_24hr', '5uM_48hr']
features = ['log2fc_hoechst_mean', 'delta_sum_mean']

# PARAMETERS
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
filters = ['cc_filter']
filters1 = []

df = pd.DataFrame()
for batch in batches:
    data = pd.read_csv('%s/01_summary/%s_grouping.txt' % (data_dir, batch), na_values=['.'], sep='\t').fillna(0)
    df['group'] = data['group']
    df['treatment'] = data['treatment']
    for feature in features:
        df['%s_%s' % (batch, feature)] = data[feature]
    if len(filters) != 0:
        for f in filters:
            df['%s_%s' % (batch, f)] = data[f]
    if len(filters1) != 0:
        for f in filters1:
            df['%s_%s' % (batch, f)] = data[f]

print(len(df))
df_flt = uti.get_df_flt(df, filters, filters1, batches)
print(len(df_flt))
df_t = uti.get_df_t(df_flt, features, batches)

fig, ax = plt.subplots(figsize=(9, 7))
# fig.subplots_adjust(right=0.8)
sns.scatterplot(data=df_t, x=features[0], y=features[1], s=8, alpha=1, color='#cccccc')  ##4d4e4e
sns.scatterplot(data=df_t[df_t['treatment'] == 'DMSO'], x=features[0], y=features[1], s=8, alpha=1, color=rainboo_colors[6])
# plt.savefig('%s/cell-cycle_vs_survival_update1.pdf' % (output_dir))
plt.show()

corr, _ = pearsonr(df_t[features[0]], df_t[features[1]])
print('Pearsons correlation: %.3f' % corr)