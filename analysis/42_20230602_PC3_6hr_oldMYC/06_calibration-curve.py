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
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230602_analysis_PC3DMandHSR_1uM_6hr_oldMYC/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

plate = 'PC3HSR_6hr'
cal = pd.read_csv('%s%s_calibration.txt' % (data_dir, plate), na_values=['.'], sep='\t')
seq = pd.read_csv('%s%s_sequence.txt' % (data_dir, plate), na_values=['.'], sep='\t')

cal_curve = pd.DataFrame()

t_total = max(cal['time'])
nuclear_b_lst = []
nuclear_a_lst = []
MYC_b_lst = []
MYC_a_lst = []
start_seq_lst = []
end_seq_lst = []
for i in range(t_total-1):
    cal_t1 = cal[cal['time'].isin([i+1])].copy().reset_index(drop=True)
    cal_t2 = cal[cal['time'].isin([i+2])].copy().reset_index(drop=True)
    cal_t1_nuclear = np.array(cal_t1['bg_hoechst_mean'].tolist() + cal_t1['mean_int_nuclear'].tolist())
    cal_t1_MYC = np.array(cal_t1['bg_MYC_mean'].tolist() + cal_t1['mean_int_MYC'].tolist())
    cal_t2_nuclear = np.array(cal_t2['bg_hoechst_mean'].tolist() + cal_t2['mean_int_nuclear'].tolist())
    cal_t2_MYC = np.array(cal_t2['bg_MYC_mean'].tolist() + cal_t2['mean_int_MYC'].tolist())
    model = LinearRegression().fit(cal_t1_nuclear.reshape((-1, 1)), cal_t2_nuclear)
    nuclear_b_lst.append(model.intercept_)
    nuclear_a_lst.append(model.coef_[0])
    model = LinearRegression().fit(cal_t1_MYC.reshape((-1, 1)), cal_t2_MYC)
    MYC_b_lst.append(model.intercept_)
    MYC_a_lst.append(model.coef_[0])
    start_seq_lst.append(np.mean(cal_t1['seq']))
    end_seq_lst.append(np.mean(cal_t2['seq']))

cal_curve['nuclear_b'] = nuclear_b_lst
cal_curve['nuclear_a'] = nuclear_a_lst
cal_curve['MYC_b'] = MYC_b_lst
cal_curve['MYC_a'] = MYC_a_lst
cal_curve['start_seq'] = start_seq_lst
cal_curve['end_seq'] = end_seq_lst

cal_curve.to_csv('%s%s_cal_ab.txt' % (output_dir, plate), index=False, sep='\t')
print("DONE!")




