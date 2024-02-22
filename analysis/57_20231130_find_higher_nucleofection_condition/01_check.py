import pandas as pd
import shared.dataframe as dat
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import seaborn as sns
import napari
import nd2

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231130_DM_nucleofection/"

data = pd.read_csv('%s/summary_treatment.txt' % master_folder, na_values=['.'], sep='\t')
print(len(data))
data['intensity'] = [dat.str_to_float(data['intensity'][i]) for i in range(len(data))]
data['positive'] = [dat.str_to_float(data['positive'][i]) for i in range(len(data))]
data['positive_average'] = [np.mean(data['positive'][i]) if len(data['positive'][i]) > 0 else 0 for i in range(len(data))]
data['positive_std'] = [np.std(data['positive'][i]) if len(data['positive'][i]) > 0 else 0 for i in range(len(data))]

print(len(data))

data.to_csv('%s/summary_treatment_new.txt' % master_folder, index=False, sep='\t')