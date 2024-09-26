

import pandas as pd
import os

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
batches = ['point5uM_24hr', 'point5uM_48hr', '5uM_24hr', '5uM_48hr']
feature = 'cytokinesis_hit'

# PARAMETERS
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']

df = pd.DataFrame()
for batch in batches:
    data = pd.read_csv('%s/01_summary/%s_grouping.txt' % (data_dir, batch), na_values=['.'], sep='\t')
    df['group'] = data['group']
    df['treatment'] = data['treatment']
    df['%s_hit' % batch] = data[feature]
    df['%s_filter' % batch] = data['qc_filter']

hit_list = []
for i in range(len(df)):
    sum_filter = 0
    sum_hit = 0
    for batch in batches:
        sum_filter = sum_filter + df['%s_filter' % batch][i]
        sum_hit = sum_hit + df['%s_hit' % batch][i]
    if sum_filter != len(batches):
        hit_list.append(-1)
    elif sum_hit >= 1:
        hit_list.append(1)
    else:
        hit_list.append(0)

df_hit = pd.DataFrame()
df_hit['group'] = df['group']
df_hit['treatment'] = df['treatment']
df_hit['hit'] = hit_list
if not os.path.exists("%s/03_hit/" % output_dir):
    os.makedirs("%s/03_hit/" % output_dir)
df_hit.to_csv('%s/03_hit/02_cytokinesis.txt' % (output_dir), index=False, sep='\t')
print("DONE!")
