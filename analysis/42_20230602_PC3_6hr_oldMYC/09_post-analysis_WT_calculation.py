import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
from shared.sinaplot import sinaplot
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230602_analysis_PC3DMandHSR_1uM_6hr_oldMYC/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_nuclear_convex_dilation = -3

# samples
exp = 'PC3HSR_6hr'
ctrl = pd.read_csv('%s%s_ctrl.txt' % (data_dir, exp), na_values=['.'], sep='\t')

# load data
data_WT = pd.DataFrame()
for i in range(len(ctrl)):
    plate = ctrl['plate'][i]
    df = pd.read_csv('%s%s_n%s.txt' % (data_dir1, plate, n_nuclear_convex_dilation), na_values=['.'], sep='\t')
    sample = ctrl['sample'][i]
    print("%s: %s" % (plate, sample))
    df_temp = df[df['sample'] == sample].copy().reset_index(drop=True)
    data_WT = pd.concat([data_WT, df_temp], axis=0)
data_WT = data_WT.reset_index(drop=True)
data_WT.to_csv('%s%s_n%s_WT.txt' % (output_dir, exp, n_nuclear_convex_dilation), index=False, sep='\t')

print("DONE!")

