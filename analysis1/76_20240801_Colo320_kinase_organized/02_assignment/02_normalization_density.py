import numpy as np
import pandas as pd

# INPUT PARAMETERS
# file info
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
data_dir1 = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/notes/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
batch = '24hr_density'

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

i = 0
start_sample = acquire_batch['sample'][i]
start_sample_index = samples_lst.index(start_sample)
n_sample = acquire_batch['n'][i]
samples_acquire = samples_lst[start_sample_index:(start_sample_index + n_sample)]
plate_acquire = acquire_batch['plate'][i]
data_acquire = data[(data['plate'] == plate_acquire) & (data['sample'].isin(samples_acquire))].copy().reset_index(drop=True)
pd_skip_acquire = pd_skip_samples[(pd_skip_samples['batch'] == batch) & (pd_skip_samples['plate'] == plate_acquire) & (pd_skip_samples['sample'].isin(samples_acquire))].copy().reset_index(drop=True)
if len(pd_skip_acquire) != 0:
    data_acquire = data_acquire[~data_acquire['sample'].isin(pd_skip_acquire['sample'].tolist())].copy().reset_index(drop=True)

data_acquire['log10_fov_hoechst_normalized'] = data_acquire['log10_fov_hoechst']
data_acquire['n_filtered_normalized'] = data_acquire['n_filtered']
data_acquire['n_pos_normalized'] = data_acquire['n_pos']
data_acquire['n_neg_normalized'] = data_acquire['n_neg']
data_acquire['acquire'] = [i+1] * len(data_acquire)

data_acquire.to_csv('%s/01_summary/%s_normalize.txt' % (output_dir, batch), index=False, sep='\t')

print("DONE!")








