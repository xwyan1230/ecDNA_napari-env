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

batch = '48hr_density'
samples =['XY%s' % (x+201) for x in range(10)] + ['XY%s' % (x+121) for x in range(10)]
identities = ['neg', 'pos'] * 5 + ['mix'] * 10
# colors = ['#808080', '#bc4d4a'] * 5 + []

df = pd.DataFrame(columns=['sample', 'identity', 'n_filtered', 'percentage'])
for i in range(len(samples)):
    sample = samples[i]
    data = pd.read_csv('%s/%s/%s/txt/%s.txt' % (data_dir, batch, batch, sample), na_values=['.'], sep='\t')
    data['log10_emiRFP670'] = np.log10(data['emiRFP670'])
    data_filtered = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])]
    n_filtered = len(data_filtered)
    n_neg = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] < cutoff)])
    percent = n_neg/n_filtered
    print(percent)
    df.loc[len(df.index)] = [sample, identities[i], n_filtered, percent]
    plt.subplots(figsize=(9, 7))
    plt.hist(data_filtered['log10_emiRFP670'], weights=np.ones(len(data_filtered)) / len(data_filtered), range=[2.5, 5], bins=40, color='w', edgecolor='black')
    plt.axvline(x=cutoff, color='#bc4d4a', linestyle='--')
    # ['#E7EAF7',(0.85,0.35,0.25)]
    # '#7BA0FB'
    plt.xlim([2.5, 5])
    plt.ylim([0, 0.6])
    if not os.path.exists('%s/%s/' % (output_dir, batch)):
        os.makedirs('%s/%s/' % (output_dir, batch))
    plt.savefig('%s/%s/%s_cutoff_qc_%s.pdf' % (output_dir, batch, batch, sample))
    plt.show()

df.to_csv('%s/%s/%s_cutoff_qc.txt' % (data_dir, batch, batch), index=False, sep='\t')

"""plt.subplots(figsize=(9, 9))
sns.violinplot(data=df[df['n_filtered'] > 500], x='identity', y='percentage')
sns.swarmplot(data=df[df['n_filtered'] > 500], x='identity', y='percentage', color="white")
plt.ylim([0, 1.05])
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_cutoff_qc.pdf' % (output_dir, batch, batch))
plt.show()"""

