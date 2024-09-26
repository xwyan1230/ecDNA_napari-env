import pandas as pd
import os

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
batches = ['point5uM_24hr', 'point5uM_48hr', '5uM_24hr', '5uM_48hr']

df = pd.DataFrame()
for batch in batches:
    data = pd.read_csv('%s/01_summary/%s_grouping.txt' % (data_dir, batch), na_values=['.'], sep='\t')
    df['group'] = data['group']
    df['treatment'] = data['treatment']
    df['%s_filter' % batch] = data['qc_filter']

df_group_seq = pd.read_csv('%s/02_group_seq/01_survival_log2fc_hoechst.txt' % data_dir, na_values=['.'], sep='\t')
hit_group_1 = df_group_seq['group'].tolist()[:20]
hit_group_2 = df_group_seq['group'].tolist()[20:37]

hit_list = []
for i in range(len(df)):
    sum_filter = 0
    for batch in batches:
        sum_filter = sum_filter + df['%s_filter' % batch][i]
    if sum_filter != len(batches):
        hit_list.append(-1)
    elif df['group'][i] in hit_group_1:
        hit_list.append(1)
    elif df['group'][i] in hit_group_2:
        hit_list.append(2)
    else:
        hit_list.append(0)

df_hit = pd.DataFrame()
df_hit['group'] = df['group']
df_hit['treatment'] = df['treatment']
df_hit['hit'] = hit_list
if not os.path.exists("%s/03_hit/" % output_dir):
    os.makedirs("%s/03_hit/" % output_dir)
df_hit.to_csv('%s/03_hit/01_survival_hoechst.txt' % (output_dir), index=False, sep='\t')
print("DONE!")



