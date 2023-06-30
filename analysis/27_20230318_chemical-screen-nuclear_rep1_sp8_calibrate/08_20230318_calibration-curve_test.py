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
output_dir = "%sdata/" % master_folder

plate = 'HSR_6hr'
cal = pd.read_csv('%s%s_calibration.txt' % (data_dir, plate), na_values=['.'], sep='\t')
cal_curve = pd.read_csv('%s%s_cal_ab.txt' % (data_dir, plate), na_values=['.'], sep='\t')

test_ch = 'MYC'
test_int = 20000
test_seq = 30
test_range = dat.find_pos(test_seq, cal_curve['start_seq'].tolist())
if test_range > 0:
    t1 = cal_curve['start_seq'][test_range-1]
    t2 = cal_curve['end_seq'][test_range-1]
    a = cal_curve['%s_a' % test_ch][test_range-1]
    b = cal_curve['%s_b' % test_ch][test_range-1]
    test_int_cal_to_region = (b * (test_seq - t1) + test_int * (t1 - t2))/(a * (t1 - test_seq) + test_seq - t2)
    test_int_initial = test_int_cal_to_region
    j = test_range
    while j > 1:
        test_int_initial = (test_int_initial - cal_curve['%s_b' % test_ch][j-2])/cal_curve['%s_a' % test_ch][j-2]
        j = j-1

plt.subplots(figsize=(9, 6))
t_total = max(cal['time'])
for i in range(t_total):
    cal_t1 = cal[cal['time'].isin([i+1])].copy().reset_index(drop=True)
    cal_t1_MYC = cal_t1['bg_MYC_mean'].tolist() + cal_t1['mean_int_MYC'].tolist()
    cal_t1_seq = cal_t1['seq'].tolist() + cal_t1['seq'].tolist()
    plt.scatter(cal_t1_seq, cal_t1_MYC, color='black')

plt.scatter([1, test_seq], [test_int_initial, test_int], color='red')
plt.show()