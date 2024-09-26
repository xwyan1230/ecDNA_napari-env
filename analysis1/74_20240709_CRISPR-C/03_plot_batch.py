import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import skimage.io as skio
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240709_analysis_CRISPR-C/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

hc = [5000, 40000]
cutoff = 3.65

samples =['XY0%s' % (x+1) for x in range(8)]
identities = ['sg2'] * 4 + ['ctrl'] * 4

df = pd.DataFrame(columns=['sample', 'identity', 'n_filtered', 'percentage'])
for i in range(len(samples)):
    sample = samples[i]
    data = pd.read_csv('%s/txt/%s.txt' % (data_dir, sample), na_values=['.'], sep='\t')
    data['log10_GFP'] = np.log10(data['GFP'])
    data_filtered = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])]
    n_filtered = len(data_filtered)
    n_pos = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_GFP'] > cutoff)])
    percent = n_pos/n_filtered
    print(percent)
    df.loc[len(df.index)] = [sample, identities[i], n_filtered, percent]
    plt.subplots(figsize=(9, 7))
    plt.hist(data_filtered['log10_GFP'], weights=np.ones(len(data_filtered)) / len(data_filtered), range=[2.5, 5],
             bins=40, color='w', edgecolor='black')
    plt.axvline(x=cutoff, color='#bc4d4a', linestyle='--')
    plt.xlim([3.0, 4.6])
    plt.ylim([0, 0.5])
    if not os.path.exists('%s/%s/' % (output_dir, sample)):
        os.makedirs('%s/%s/' % (output_dir, sample))
    plt.savefig('%s/%s/%s_cutoff_qc.pdf' % (output_dir, sample, sample))
    plt.show()

df.to_csv('%s/cutoff_qc.txt' % (data_dir), index=False, sep='\t')

plt.subplots(figsize=(9, 7))
sns.barplot(data=df[df['n_filtered'] > 200], x='identity', y='percentage', color='#e3e3e3', edgecolor='#4d4e4e')
sns.swarmplot(data=df[df['n_filtered'] > 200], x='identity', y='percentage', color='#bc4d4a')
plt.ylim([0, 0.05])
plt.savefig('%s/cutoff_qc.pdf' % output_dir)
plt.show()

