import pandas as pd
import shared.display as dis
import numpy as np
import os
import matplotlib.pyplot as plt

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220809_Aarohi_chromosomeRecombination_metaphase" \
                "/08052022_ChrRecomb_KCl_timept/"
sample = 'ATCC'
save_folder = master_folder

data = pd.read_csv(("%s%s.txt" % (master_folder, sample)), na_values=['.'], sep='\t')

print(len(data))

plt.scatter(data['DM_copy'], data['HSR_copy'], s=5)
plt.xlabel('DM copy number')
plt.ylabel('HSR copy number')
plt.show()