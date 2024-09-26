import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.stats import pearsonr
import skimage.io as skio
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240227_analysis_SNU16_co-segregation/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

data = pd.read_csv('%s/Cosegregation_Summary.csv' % data_dir, na_values=['.'], sep=',')
data_ctrl = data[data['Treatment'] == 'Control'].copy().reset_index(drop=True)
data_wee1 = data[data['Treatment'] == 'dWEE1'].copy().reset_index(drop=True)

n1 = len(data_ctrl)
n2 = len(data_wee1)
print(n1)
print(n2)

# calculate pearson correlation R
r1, p1 = pearsonr(data_ctrl['Percent_Gene_1_Inherited'], data_ctrl['Percent_Gene_2_Inherited'])
print('Pearsons correlation: %.3f, %s' % (r1, p1))
r2, p2 = pearsonr(data_wee1['Percent_Gene_1_Inherited'], data_wee1['Percent_Gene_2_Inherited'])
print('Pearsons correlation: %.3f, %s' % (r2, p2))

# Fisher's z test
r1_trans = 0.5 * np.log((1+r1)/(1-r1))
r2_trans = 0.5 * np.log((1+r2)/(1-r2))
z = (r1_trans - r2_trans)/(np.sqrt((1/(n1-3))+(1/(n2-3))))

print(r1_trans)
print(r2_trans)
print(z)
# https://www.socscistatistics.com/pvalues/normaldistribution.aspx
p = 0.055827  # single sided

