import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

# INPUT PARAMETERS
# file info
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
batch = '5uM_24hr'

# PARAMETERS
data = pd.read_csv('%s/01_summary/%s.txt' % (data_dir, batch), na_values=['.'], sep='\t')

fig, ax = plt.subplots(figsize=(9, 7))
# fig.subplots_adjust(right=0.8)
sns.scatterplot(data=data, x='n_filtered', y='check_green', s=8, alpha=1, color='#cccccc')
print(data[data['check_green'] > 50]['plate'])
print(data[data['check_green'] > 50]['sample'])
print(data[data['check_green'] > 50]['treatment'])
plt.show()