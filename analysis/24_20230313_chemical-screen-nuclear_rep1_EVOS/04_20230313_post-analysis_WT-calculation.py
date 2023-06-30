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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230313_analysis_chemical-screen-nuclear_rep1_EVOS/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

# samples
exp = 'HSR_6hr'
ctrls = ['HSR_6hr*C3', 'HSR_6hr*C10', 'HSR_6hr*D6', 'HSR_6hr*F3', 'HSR_6hr*F10']

# load data
data_WT = pd.DataFrame()
for i in range(len(ctrls)):
    plate = ctrls[i].split('*')[0]
    df = pd.read_csv('%s%s.txt' % (data_dir, plate), na_values=['.'], sep='\t')
    well = ctrls[i].split('*')[1]
    print("%s: %s" % (plate, well))
    df_temp = df[df['well'] == well].copy().reset_index(drop=True)
    data_WT = pd.concat([data_WT, df_temp], axis=0)
data_WT = data_WT.reset_index(drop=True)
data_WT.to_csv('%s%s_WT.txt' % (output_dir, exp), index=False, sep='\t')

print("DONE!")

