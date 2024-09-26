# Fig 2e

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
batches = ['point5uM_24hr', 'point5uM_48hr', '5uM_24hr', '5uM_48hr']

# PARAMETERS
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
mode = 'neg'  # only accepts 'neg' and 'pos'
feature = 'log2fc_n_%s_mean' % mode

if mode == 'neg':
    color_theme = [9, 8, 7, 6, 5]
else:
    color_theme = [0, 1, 2, 3, 4]

df = pd.DataFrame()
for batch in batches:
    data = pd.read_csv('%s/01_summary/%s_hitcalling.txt' % (data_dir, batch), na_values=['.'], sep='\t')
    data1 = pd.read_csv('%s/01_summary/%s_grouping.txt' % (data_dir, batch), na_values=['.'], sep='\t')
    df['group'] = data['group']
    df['treatment'] = data['treatment']
    df['%s_%s_hit' % (batch, feature)] = data[feature]
    df['%s_%s' % (batch, feature)] = data1[feature]


for batch in batches:
    df = df[df['%s_%s_hit' % (batch, feature)] != -100].copy().reset_index(drop=True)
print(len(df))

df['hit_sum'] = np.sum(df[['%s_%s_hit' % (batch, feature) for batch in batches]], axis=1)
df['sum'] = np.sum(df[['%s_%s' % (batch, feature) for batch in batches]], axis=1)

df_num = df.copy()
df_num = df_num[df_num['hit_sum'] <= -1]
df_num = df_num.sort_values(by=['hit_sum', 'sum'], ascending=[True, True]).reset_index(drop=True)

n_total = len(df[df['treatment'] != 'DMSO'])
print(n_total)
n4_total = len(df[df['hit_sum'] == -4])
n4 = len(df[(df['hit_sum'] == -4) & (df['treatment'] != 'DMSO')])
n3_total = len(df[df['hit_sum'] == -3])
n3 = len(df[(df['hit_sum'] == -3) & (df['treatment'] != 'DMSO')])
n2_total = len(df[df['hit_sum'] == -2])
n2 = len(df[(df['hit_sum'] == -2) & (df['treatment'] != 'DMSO')])
n1_total = len(df[df['hit_sum'] == -1])
n1 = len(df[(df['hit_sum'] == -1) & (df['treatment'] != 'DMSO')])
n0 = len(df[(df['hit_sum'] >= 0) & (df['treatment'] != 'DMSO')])

print("4: %s/%s" % (n4, n4_total))
print("3: %s/%s" % (n3, n3_total))
print("2: %s/%s" % (n2, n2_total))
print("1: %s/%s" % (n1, n1_total))

plt.subplots(figsize=(9, 7))
y = np.array([n4, n3, n2, n1, n0])
print(y)
plt.pie(y, colors=[rainboo_colors[color_theme[0]], rainboo_colors[color_theme[1]], rainboo_colors[color_theme[2]],
                   rainboo_colors[color_theme[3]], rainboo_colors[color_theme[4]]])
# plt.savefig('%s/pie_hoechst_neg_hit.pdf' % output_dir)
plt.show()

if not os.path.exists("%s/04_others/" % output_dir):
    os.makedirs("%s/04_others/" % output_dir)
df.to_csv('%s/04_others/01_survival_n_%s.txt' % (output_dir, mode), index=False, sep='\t')