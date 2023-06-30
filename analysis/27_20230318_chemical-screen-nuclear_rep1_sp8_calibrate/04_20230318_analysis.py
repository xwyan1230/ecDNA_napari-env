import skimage.io as skio
import numpy as np
from skimage.measure import regionprops
from skimage.morphology import erosion, disk
import math
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230318_analysis_chemical-screen-nuclear_rep1_sp8_calibrate/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_nuclear_convex_dilation = -3

data = pd.DataFrame(columns=['plate', 'sample', 'well', 'nuclear', 'FOV', 'n_nuclear_convex_dilation',
                             'centroid_nuclear', 'area_nuclear', 'circ_nuclear', 'mean_int_nuclear', 'mean_int_MYC'])

plate = 'DM_6hr'
seq = pd.read_csv('%s%s_sequence.txt' % (data_dir, plate), na_values=['.'], sep='\t')

samples = seq['sample'].tolist()
start_sample = 0

for s in range(len(samples)):
    print(s+start_sample)
    sample = samples[s+start_sample]
    print(sample)

    # file_name = '20230318_HSR_6hr_rep1_calibration_%s_RAW' % (sample)
    file_name = '20230319_DM_6hr_rep1_calibrate_%s_RAW' % (sample)
    img_hoechst_stack = skio.imread("%s%s/%s_ch01.tif" % (data_dir, plate, file_name), plugin="tifffile")
    img_MYC_stack = skio.imread("%s%s/%s_ch00.tif" % (data_dir, plate, file_name), plugin="tifffile")

    for fov in range(4):
        img_hoechst = img_hoechst_stack[:, :, fov]
        img_MYC = img_MYC_stack[:, :, fov]
        # img_MYC = np.concatenate([np.zeros(shape=[5, 2048]), img_MYC], axis=0)[:2048, :2048]
        img_nuclear_seg = skio.imread("%s%s/seg_tif/%s_seg_%s.tif" % (data_dir1, plate, file_name, fov), plugin="tifffile")
        if n_nuclear_convex_dilation < 0:
            img_nuclear_seg = erosion(img_nuclear_seg, disk(abs(n_nuclear_convex_dilation)))

        nuclear_props = regionprops(img_nuclear_seg, img_hoechst)
        MYC_props = regionprops(img_nuclear_seg, img_MYC)

        centroid_nuclear = [nuclear_props[x].centroid for x in range(len(nuclear_props))]
        area_nuclear = [nuclear_props[x].area for x in range(len(nuclear_props))]
        perimeter_nuclear = [nuclear_props[x].perimeter for x in range(len(nuclear_props))]
        mean_int_nuclear = [nuclear_props[x].intensity_mean for x in range(len(nuclear_props))]
        mean_int_IF = [MYC_props[x].intensity_mean for x in range(len(MYC_props))]
        circ_nuclear = [(4 * math.pi * area_nuclear[x]) / (perimeter_nuclear[x] ** 2) for x in range(len(area_nuclear))]

        data_temp = pd.DataFrame({'plate': [plate] * len(area_nuclear),
                                  'sample': [sample] * len(area_nuclear),
                                  'well': [sample.split('_')[0]] * len(area_nuclear),
                                  'nuclear': range(len(area_nuclear)),
                                  'FOV': [fov] * len(area_nuclear),
                                  'n_nuclear_convex_dilation': [-3] * len(area_nuclear),
                                  'centroid_nuclear': centroid_nuclear,
                                  'area_nuclear': area_nuclear,
                                  'circ_nuclear': circ_nuclear,
                                  'mean_int_nuclear': mean_int_nuclear,
                                  'mean_int_MYC': mean_int_IF})
        data = pd.concat([data, data_temp], axis=0)

data.to_csv('%s%s_n%s.txt' % (output_dir, plate, n_nuclear_convex_dilation), index=False, sep='\t')
print("DONE!")
