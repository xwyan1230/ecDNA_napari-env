import pandas as pd
import shared.display as dis
import numpy as np
import os
import matplotlib.pyplot as plt

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220809_Aarohi_chromosomeRecombination_metaphase" \
                "/08052022_ChrRecomb_KCl_timept/"
save_folder = master_folder

sample = 'Jun16'
get_rid_fov = [37, 46, 52]

sample1 = 'ATCC'
get_rid_fov1 = []

sample2 = 'XY'
get_rid_fov2 = []

data = pd.read_csv(("%s%s.txt" % (master_folder, sample)), na_values=['.'], sep='\t')
data_filter = data[~data['FOV'].isin(get_rid_fov)]
data1 = pd.read_csv(("%s%s.txt" % (master_folder, sample1)), na_values=['.'], sep='\t')
data_filter1 = data1[~data1['FOV'].isin(get_rid_fov1)]
data2 = pd.read_csv(("%s%s.txt" % (master_folder, sample2)), na_values=['.'], sep='\t')
data_filter2 = data2[~data2['FOV'].isin(get_rid_fov2)]

plt.scatter(data_filter['DM_copy'], data_filter['HSR_copy'], s=5, label=sample)
plt.scatter(data_filter1['DM_copy'], data_filter1['HSR_copy'], s=5, label=sample1)
plt.scatter(data_filter2['DM_copy'], data_filter2['HSR_copy'], s=5, label=sample2)
plt.xlabel('DM copy number')
plt.ylabel('HSR copy number')
plt.legend()
plt.ylim([0, 500])
plt.show()