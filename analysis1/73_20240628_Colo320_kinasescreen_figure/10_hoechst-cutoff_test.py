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
cutoff = 3.1

folder = '5uM_48hr'
batch = '5uM_48hr_1'
if int(batch[-1]) == 1:
    samples = ['XY02', 'XY04', 'XY06', 'XY08', 'XY10', 'XY12', 'XY14', 'XY16', 'XY18', 'XY20']
    samples = ['XY0%s' % (x+1) for x in range(9)] + ['XY%s' % (x+10) for x in range(191)]
elif int(batch[-1]) == 2:
    samples = ['XY%s' % (2*x+111) for x in range(5)] + ['XY%s' % (2*x+101) for x in range(5)]
else:
    samples = ['XY%s' % (2*x+201) for x in range(5)] + ['XY%s' % (2*x+202) for x in range(5)]
"""samples = ['XY02', 'XY04', 'XY06', 'XY08', 'XY10', 'XY12', 'XY14', 'XY16', 'XY18', 'XY20'] + \
          ['XY01', 'XY03', 'XY05', 'XY07', 'XY09', 'XY11', 'XY13', 'XY15', 'XY17', 'XY19'] + \
          ['XY%s' % (2*x+41) for x in range(10)] + ['XY%s' % (2*x+101) for x in range(10)] + \
          ['XY%s' % (2*x+182) for x in range(10)]"""
"""samples = ['XY%s' % (2*x+111) for x in range(5)] + ['XY%s' % (2*x+101) for x in range(5)] + \
          ['XY%s' % (2*x+41) for x in range(10)] + \
          ['XY%s' % (2*x+102) for x in range(10)] + ['XY%s' % (2*x+142) for x in range(10)]"""
"""samples = ['XY%s' % (2*x+201) for x in range(5)] + ['XY%s' % (2*x+202) for x in range(5)] + \
          ['XY02', 'XY04', 'XY06', 'XY08', 'XY10', 'XY12', 'XY14', 'XY16', 'XY18', 'XY20'] + \
          ['XY%s' % (2*x+81) for x in range(10)] + ['XY%s' % (2*x+142) for x in range(10)]"""

identities = ['neg'] * 5 + ['pos'] * 5 + ['mix'] * 40
# colors = ['#808080', '#bc4d4a'] * 5 + []

df = pd.DataFrame(columns=['sample', 'identity', 'n_filtered', 'percentage'])
for i in range(len(samples)):
    sample = samples[i]
    print(sample)
    data = pd.read_csv('%s/%s/%s/txt/%s.txt' % (data_dir, folder, batch, sample), na_values=['.'], sep='\t')
    data['log10_emiRFP670'] = np.log10(data['emiRFP670'])
    data_filtered = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])]
    n_filtered = len(data_filtered)
    print(n_filtered)
    n_neg = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] < cutoff)])
    percent = n_neg/n_filtered
    print(percent)
    # df.loc[len(df.index)] = [sample, identities[i], n_filtered, percent]
    plt.subplots(figsize=(9, 7))
    m = plt.hist(data_filtered['log10_emiRFP670'], weights=np.ones(len(data_filtered)) / len(data_filtered), range=[2.5, 5], bins=60, color='w', edgecolor='black')
    # print(m)
    test_index = m[0].tolist().index(np.max(m[0][:15]))
    test_val = m[1][test_index]
    print(test_index)
    print(test_val)
    plt.axvline(x=test_val+0.3, color='#bc4d4a', linestyle='--')
    # ['#E7EAF7',(0.85,0.35,0.25)]
    # '#7BA0FB'
    plt.xlim([2.5, 5])
    plt.ylim([0, 0.6])
    """if not os.path.exists('%s/%s/%s/cutoff_qc/' % (output_dir, folder, batch)):
        os.makedirs('%s/%s/%s/cutoff_qc/' % (output_dir, folder, batch))
    plt.savefig('%s/%s/%s/cutoff_qc/%s_cutoff_qc_%s.pdf' % (output_dir, folder, batch, batch, sample))"""
    plt.show()
    # plt.close()

# df.to_csv('%s/%s/%s/%s_cutoff_qc.txt' % (data_dir, folder, batch, batch), index=False, sep='\t')

"""plt.subplots(figsize=(9, 9))
sns.violinplot(data=df[df['n_filtered'] > 500], x='identity', y='percentage')
sns.swarmplot(data=df[df['n_filtered'] > 500], x='identity', y='percentage', color="white")
plt.ylim([0, 1.05])
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_cutoff_qc.pdf' % (output_dir, batch, batch))
plt.show()"""

