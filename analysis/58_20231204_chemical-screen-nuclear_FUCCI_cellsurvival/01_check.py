import pandas as pd
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import seaborn as sns

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231204_analysis_FUCCI-screen/"
sample = 'HSR_FUCCI_6hr'
ctrl = ['C3', 'C10', 'D6', 'F3', 'F10']

data = pd.read_csv('%s/%s/summary.txt' % (master_folder, sample), na_values=['.'], sep='\t')
data['total'] = data['G1'] + data['S/G2'] + data['G2/M'] + data['G1/S']

avg_ctrl = np.mean(data[data['sample'].isin(ctrl)]['total'].tolist())
data['per_survival'] = data['total']/avg_ctrl
data['per_survival_log2FC'] = np.log2(data['per_survival'])

data.to_csv('%s/%s/summary_survival.txt' % (master_folder, sample), index=False, sep='\t')

