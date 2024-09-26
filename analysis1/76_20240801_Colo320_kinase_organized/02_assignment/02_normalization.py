import numpy as np
import pandas as pd
import utilities as uti
import os

# INPUT PARAMETERS
# filter out:
# 1. wells that miss 30% images
# 2. DMSO wells that have less than 1000 identified cells
# file info
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
data_dir1 = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/notes/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
batch = '5uM_48hr'

# PARAMETERS
n_img = 9
skip_cutoff = 0.3
skip = pd.read_csv('%s/skip.txt' % data_dir1, na_values=['.'], sep='\t')
skip_fovs = [i.split('-')[0] + '-' + i.split('-')[1] + '-' + i.split('-')[2] for i in skip['samples'].tolist()]
skip_fovs_dict = {i: skip_fovs.count(i) for i in skip_fovs}
pd_skip_fovs = pd.DataFrame(skip_fovs_dict.items(), columns=['samples', 'n'])
pd_skip_fovs['batch'] = [i.split('-')[0] for i in pd_skip_fovs['samples'].tolist()]
pd_skip_fovs['plate'] = [i.split('-')[1] for i in pd_skip_fovs['samples'].tolist()]
pd_skip_fovs['plate'] = pd_skip_fovs['plate'].astype(int)
pd_skip_fovs['sample'] = [i.split('-')[2] for i in pd_skip_fovs['samples'].tolist()]
pd_skip_samples = pd_skip_fovs[pd_skip_fovs['n'] > (n_img * skip_cutoff)].copy().reset_index(drop=True)

data = pd.read_csv('%s/01_summary/%s.txt' % (data_dir, batch), na_values=['.'], sep='\t')
data['log10_fov_hoechst'] = np.log10(data['fov_hoechst'])
acquire = pd.read_csv('%s/acquire.txt' % data_dir1, na_values=['.'], sep='\t')
acquire_batch = acquire[acquire['batch'] == batch].copy().reset_index(drop=True)
samples_lst = ['XY0%s' % (x + 1) for x in range(9)] + ['XY%s' % (x + 10) for x in range(201)]

acquire_lst = []
data_acquire_lst = []
data_acquire_DMSO_lst = []
len_lst = []

for i in range(len(acquire_batch)):
    start_sample = acquire_batch['sample'][i]
    start_sample_index = samples_lst.index(start_sample)
    n_sample = acquire_batch['n'][i]
    samples_acquire = samples_lst[start_sample_index:(start_sample_index + n_sample)]
    plate_acquire = acquire_batch['plate'][i]
    data_acquire = data[(data['plate'] == plate_acquire) & (data['sample'].isin(samples_acquire))].copy().reset_index(drop=True)
    data_acquire_filter = data_acquire[(data_acquire['treatment'] == 'DMSO') & (data_acquire['n_filtered'] <= 1000)]
    data_acquire = data_acquire.drop(data_acquire_filter.index).reset_index(drop=True)
    pd_skip_acquire = pd_skip_samples[(pd_skip_samples['batch'] == batch) & (pd_skip_samples['plate'] == plate_acquire) & (pd_skip_samples['sample'].isin(samples_acquire))].copy().reset_index(drop=True)
    if len(pd_skip_acquire) != 0:
        data_acquire = data_acquire[~data_acquire['sample'].isin(pd_skip_acquire['sample'].tolist())].copy().reset_index(drop=True)
    data_acquire_DMSO = data_acquire[data_acquire['treatment'] == 'DMSO'].copy().reset_index(drop=True)
    if len(data_acquire_DMSO) > 5:
        data_acquire_lst.append(data_acquire)
        data_acquire_DMSO_lst.append(data_acquire_DMSO)
        acquire_lst.append(i+1)
        len_lst.append(len(data_acquire_DMSO))

normalize_aquire = len_lst.index(np.max(len_lst)) + 1
ref = data_acquire_DMSO_lst[len_lst.index(np.max(len_lst))]
data_normalize = pd.DataFrame()
for i in range(len(acquire_lst)):
    data_acquire = uti.get_normalization(data_acquire_lst[i], data_acquire_DMSO_lst[i], ref)
    data_acquire['acquire'] = [acquire_lst[i]] * len(data_acquire)
    data_normalize = pd.concat([data_normalize, data_acquire], axis=0)

data_normalize.to_csv('%s/01_summary/%s_normalize.txt' % (output_dir, batch), index=False, sep='\t')

if os.path.exists("%s/normalization.txt" % data_dir1):
    norm = pd.read_csv('%s/normalization.txt' % data_dir1, na_values=['.'], sep='\t')
else:
    norm = pd.DataFrame(columns=['batch', 'acquire'])
if len(norm[norm['batch'] == batch]) == 0:
    norm.loc[len(norm.index)] = [batch, normalize_aquire]
else:
    location = norm[norm['batch'] == batch].index
    norm.loc[location] = [batch, normalize_aquire]
norm.to_csv('%s/normalization.txt' % data_dir1, index=False, sep='\t')
print("DONE!")








