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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230602_analysis_PC3DMandHSR_1uM_6hr_oldMYC/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

plate = 'PC3HSR_6hr'
cal = pd.read_csv('%s%s_calibration.txt' % (data_dir, plate), na_values=['.'], sep='\t')
cal_curve = pd.read_csv('%s%s_cal_ab.txt' % (data_dir, plate), na_values=['.'], sep='\t')

test_chs = ['MYC', 'nuclear']
for test_ch in test_chs:
    test_int = 20000

    if test_ch == 'MYC':
        ch1 = 'MYC'
    else:
        ch1 = 'hoechst'

    plt.subplots(figsize=(9, 6))
    t_total = max(cal['time'])
    for i in range(t_total):
        cal_t1 = cal[cal['time'].isin([i+1])].copy().reset_index(drop=True)
        cal_t1_MYC = cal_t1['bg_%s_mean' % ch1].tolist() + cal_t1['mean_int_%s' % test_ch].tolist()
        plt.scatter([i+1]*len(cal_t1_MYC), cal_t1_MYC, color='black')

    plt.scatter([1], [test_int], color='red')
    for i in range(t_total-1):
        test_int = cal_curve['%s_a' % test_ch][i]*test_int+cal_curve['%s_b' % test_ch][i]
        plt.scatter([i+2], [test_int], color='red')
    plt.savefig('%s/%s_calibration_%s.pdf' % (output_dir, plate, test_ch))
    plt.show()