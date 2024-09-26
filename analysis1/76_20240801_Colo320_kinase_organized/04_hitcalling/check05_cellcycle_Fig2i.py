import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import utilities as uti

# INPUT PARAMETERS
# file info
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"

batches = ['point5uM_24hr', 'point5uM_48hr', '5uM_24hr', '5uM_48hr']
mode = 'all'  # only accepts 'neg', 'pos', 'all'

# PARAMETERS
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']

if mode == 'all':
    features = ['log2fc_hoechst_mean', 'per_G1_mean', 'per_G1S_mean', 'per_S_mean', 'per_G2M_mean', 'per_G2MG1_mean']
    filters = ['cc_filter']
    filters1 = []
    interval = [-0.175, -0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0]
else:
    features = ['log2fc_n_%s_mean' % mode, 'per_%s_G1_mean' % mode, 'per_%s_G1S_mean' % mode, 'per_%s_S_mean' % mode,
                'per_%s_G2M_mean' % mode, 'per_%s_G2MG1_mean' % mode]
    filters = ['cc_filter']
    filters1 = ['cytokinesis_hit']
    interval = [-4.9, -4.2, -3.5, -2.8, -2.1, -1.4, -0.7, 0]

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
print(len(df_t))

survival_group = []
pointer = 0
for i in range(len(df_t)):
    for j in range(len(interval)):
        if (df_t[features[0]][i] < interval[j]) & (pointer == 0):
            survival_group.append(j)
            pointer = 1
    if pointer == 0:
        survival_group.append(8)
    else:
        pointer = 0
df_t['survival_group'] = survival_group
df_t.to_csv('%s/temp.txt' % output_dir, index=False, sep='\t')

cc = features[1:6]
df_cc = pd.DataFrame()
df_cc['survival_group'] = np.arange(0, 9, 1)
for j in range(len(cc)):
    avg_lst = []
    for i in range(9):
        avg_lst.append(np.mean(df_t[df_t['survival_group'] == i][cc[j]]))
    df_cc[cc[j]] = avg_lst
df_cc.columns = ['survival_group', 'G1', 'G1S', 'S', 'G2M', 'G2MG1']

colors = ['#bc4d4a', '#3d5a80', '#669daa', '#dd933c', '#d2b48c'] ## Colors for the bars
df_cc.set_index('survival_group').plot(kind='bar', stacked=True, color=colors) ## Plot
plt.ticklabel_format(style='plain', useOffset=False, axis='y') ## No offset
plt.gca().set_ylabel("Percentage (%)")
# plt.savefig('%s/cell-cycle_group_by_survival.pdf' % (output_dir))
plt.show()