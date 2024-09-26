# Fig 1d

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
data_dir1 = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/processed/"
# output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/figures/"
batch = '5uM_48hr'
plate = 2
sample = 'XY77'

data_batch = pd.read_csv('%s/01_summary/%s.txt' % (data_dir, batch), na_values=['.'], sep='\t')
temp = data_batch[(data_batch['plate'] == plate) & (data_batch['sample'] == sample)]
hc = [temp['hoechst_cutoff_low'].tolist()[0], temp['hoechst_cutoff_high'].tolist()[0]]
cutoff = temp['log10_emiRFP670_cutoff'].tolist()[0]
data = pd.read_csv('%s/%s/%s_%s/txt/%s.txt' % (data_dir1, batch, batch, plate, sample), na_values=['.'], sep='\t')
data['log10_emiRFP670'] = np.log10(data['emiRFP670'])
data_filtered = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])].copy().reset_index(drop=True)

plt.subplots(figsize=(9, 7))
plt.hist(data_filtered['log10_emiRFP670'], weights=np.ones(len(data_filtered)) / len(data_filtered), range=[2.5, 5], bins=40, color='w', edgecolor='black')
plt.axvline(x=cutoff, color='#bc4d4a', linestyle='--')
plt.xlim([2.5, 5])
plt.ylim([0, 0.6])
plt.show()

print("DONE!")