import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import skimage.io as skio
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
data_dir1 = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

hc = [5000, 40000]
cutoff = 3.2

df = pd.DataFrame()
batches = ['24hr_density', '48hr_density']
for batch in batches:
    data = pd.read_csv('%s/%s/%s_cell-cycle_update1.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    data_temp = pd.DataFrame()
    data_temp['group'] = [batch] * len(data) * 5
    data_temp['identity'] = data['identity'].tolist() * 5
    data_temp['cellcycle'] = ['G1'] * len(data) + ['G1S'] * len(data) + ['S'] * len(data) + ['G2M'] * len(data) + ['G2MG1'] * len(data)
    data_temp['percentage'] = data['G1'].tolist() + data['G1S'].tolist() + data['S'].tolist() + data['G2M'].tolist() + data['G2MG1'].tolist()
    data_temp['x'] = ['%s_%s' % (batch[0:4], data_temp['cellcycle'][x]) for x in range(len(data_temp))]
    df = pd.concat([df, data_temp], axis=0).reset_index(drop=True)
# df = df.sort_values(by=['cellcycle', 'group'], ascending=[False, True]).reset_index(drop=True)
df = df[df['identity'] == 'mix'].copy().reset_index(drop=True)

plt.subplots(figsize=(9, 7))
sns.barplot(data=df, x='x', y='percentage', color='#e3e3e3', edgecolor='#4d4e4e')
sns.swarmplot(data=df, x='x', y='percentage', color='#bc4d4a')
plt.ylim([0, 0.5])
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/cellcycle_update1.pdf' % output_dir)
plt.show()

