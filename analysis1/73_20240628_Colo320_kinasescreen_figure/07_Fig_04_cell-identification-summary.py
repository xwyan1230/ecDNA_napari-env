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
    data = pd.read_csv('%s/%s/%s_cutoff_qc.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    data['batch'] = [batch] * len(data)
    data['group'] = ['%s_%s' % (batch[0:4], data['identity'][x]) for x in range(len(data))]
    df = pd.concat([df, data], axis=0).reset_index(drop=True)
df = df.sort_values(by=['batch', 'identity'], ascending=[True, False]).reset_index(drop=True)

plt.subplots(figsize=(9, 7))
# plt.bar(df[df['n_filtered'] > 500]['group'], df[df['n_filtered'] > 500]['percentage'], 0.5)
sns.barplot(data=df[df['n_filtered'] > 500], x='group', y='percentage', color='#e3e3e3', edgecolor='#4d4e4e')
sns.swarmplot(data=df[df['n_filtered'] > 500], x='group', y='percentage', color='#bc4d4a')
plt.ylim([0, 1.05])
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/cutoff_qc.pdf' % output_dir)
plt.show()

