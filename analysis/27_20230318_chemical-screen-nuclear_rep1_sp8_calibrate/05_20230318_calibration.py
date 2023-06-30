import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
import shared.dataframe as dat
from skimage import segmentation
import shared.math as mat
import math
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230318_analysis_chemical-screen-nuclear_rep1_sp8_calibrate/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sdata/" % master_folder

n_nuclear_convex_dilation = -3
plate = 'DM_6hr'
cal = pd.read_csv('%s%s_calibration.txt' % (data_dir, plate), na_values=['.'], sep='\t')
seq = pd.read_csv('%s%s_sequence.txt' % (data_dir, plate), na_values=['.'], sep='\t')
df = pd.read_csv('%s%s_n%s.txt' % (data_dir1, plate, n_nuclear_convex_dilation), na_values=['.'], sep='\t')

cal['seq'] = [seq[seq['sample'] == cal['sample'][i]]['seq'].tolist()[0] for i in range(len(cal))]

samples = cal['sample'].tolist()

mean_int_nuclear_lst = []
mean_int_MYC_lst = []

for sample in samples:
    print(sample)
    df_temp = df[df['sample'] == sample].copy().reset_index(drop=True)
    mean_int_nuclear = np.mean(df_temp['mean_int_nuclear'])
    mean_int_MYC = np.mean(df_temp['mean_int_MYC'])
    mean_int_nuclear_lst.append(mean_int_nuclear)
    mean_int_MYC_lst.append(mean_int_MYC)

cal['mean_int_nuclear'] = mean_int_nuclear_lst
cal['mean_int_MYC'] = mean_int_MYC_lst

cal.to_csv('%s%s_calibration.txt' % (output_dir, plate), index=False, sep='\t')
print("DONE!")