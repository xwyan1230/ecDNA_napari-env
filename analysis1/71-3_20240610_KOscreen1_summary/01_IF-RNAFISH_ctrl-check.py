import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import skimage.io as skio
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%ssummary/" % master_folder
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

mode = 'RNAFISH'
ctrls = ['J3', 'J4', 'J5', 'J6', 'J7', 'J8', 'J9', 'J10', 'J11', 'J12']
reps = [1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
identities = ['GFP', 'mCherry'] * 5

df = pd.DataFrame(columns=['sample', 'identity', 'n_GFP', 'n_mCherry', 'percentage'])

for i in range(len(ctrls)):
    sample = ctrls[i]
    identity = identities[i]
    rep = reps[i]
    n_file = 25 if mode == 'IF' else 28
    data = pd.read_csv("%s/%s/%s_%s_%s_rep%s_red_green_group.txt" % (data_dir, sample, n_file, sample, mode, rep), na_values=['.'], sep='\t')
    n_GFP = len(data[data['group'] == 'GFP'])
    n_mCherry = len(data[data['group'] == 'mCherry'])
    if identity == 'GFP':
        percentage = n_GFP/(n_GFP+n_mCherry)
    else:
        percentage = n_mCherry/(n_GFP+n_mCherry)
    df.loc[len(df.index)] = [sample, identity, n_GFP, n_mCherry, percentage]

df.to_csv('%s/01_%s_cutoff_qc.txt' % (output_dir, mode), index=False, sep='\t')

fig, ax = plt.subplots(figsize=(4, 6))
fig.subplots_adjust(left=0.4)
sns.violinplot(data=df, x='identity', y='percentage')
sns.swarmplot(data=df, x='identity', y='percentage', color=(255 / 255, 140 / 255, 0 / 255))
plt.ylim([0, 1.05])
plt.savefig('%s/01_%s_cutoff_qc.pdf' % (output_dir, mode))
plt.show()

