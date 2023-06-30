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
from sklearn.linear_model import LinearRegression
import shared.math as mat
import math
import matplotlib.pyplot as plt
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230318_analysis_chemical-screen-nuclear_rep1_sp8_calibrate/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_nuclear_convex_dilation = -3
plate = 'DM_6hr'
cal = pd.read_csv('%s%s_calibration.txt' % (data_dir, plate), na_values=['.'], sep='\t')
cal_curve = pd.read_csv('%s%s_cal_ab.txt' % (data_dir, plate), na_values=['.'], sep='\t')
seq = pd.read_csv('%s%s_sequence.txt' % (data_dir, plate), na_values=['.'], sep='\t')
df = pd.read_csv('%s%s_n%s.txt' % (data_dir1, plate, n_nuclear_convex_dilation), na_values=['.'], sep='\t')
chs = ['nuclear', 'MYC']

calibrated_lst = []

for ch in chs:
    lst_temp = []
    for i in range(len(df)):
        test_int = df['mean_int_%s' % ch].tolist()[i]
        test_seq = seq[seq['sample'] == df['sample'][i]]['seq'].tolist()[0]
        test_range = dat.find_pos(test_seq, cal_curve['start_seq'].tolist())
        if test_range > 0:
            t1 = cal_curve['start_seq'][test_range-1]
            t2 = cal_curve['end_seq'][test_range-1]
            a = cal_curve['%s_a' % ch][test_range - 1]
            b = cal_curve['%s_b' % ch][test_range - 1]
            test_int_cal_to_region = (b * (test_seq - t1) + test_int * (t1 - t2))/(a * (t1 - test_seq) + test_seq - t2)
            test_int_initial = test_int_cal_to_region
            j = test_range
            while j > 1:
                test_int_initial = (test_int_initial - cal_curve['%s_b' % ch][j - 2]) / cal_curve['%s_a' % ch][j - 2]
                j = j-1
            lst_temp.append(test_int_initial)
        else:
            lst_temp.append(test_int)
    calibrated_lst.append(lst_temp)

for i in range(len(chs)):
    df['mean_int_%s_cal' % chs[i]] = calibrated_lst[i]

df.to_csv('%s%s_n%s.txt' % (output_dir, plate, n_nuclear_convex_dilation), index=False, sep='\t')
print("DONE!")



