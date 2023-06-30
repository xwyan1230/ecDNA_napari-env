import skimage.io as skio
import numpy as np
from skimage.measure import regionprops
import math
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230317_analysis_chemical-screen-nuclear_rep1_Mischel/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

data = pd.DataFrame(columns=['plate', 'well', 'nuclear', 'FOV', 'n_nuclear_convex_dilation',
                             'centroid_nuclear', 'area_nuclear', 'circ_nuclear', 'mean_int_nuclear', 'mean_int_MYC'])

plate = 'HSR_6hr'
rows = ['B', 'C', 'D', 'E', 'F', 'G']
columns = ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
num_total_samples = len(rows) * len(columns)
for i in range(num_total_samples):
    row = rows[int(i/len(columns))]
    column = columns[int(i - int(i/len(columns))*len(columns))]
    print('%s%s' % (row, column))
    for j in range(5):
        file_name = '%s_%s_R%s_RAW' % (row, column, j + 1)
        img_hoechst = skio.imread("%s%s/%s_ch00.tif" % (data_dir, plate, file_name), plugin="tifffile")  # 1536, 2048
        img_MYC = skio.imread("%s%s/%s_ch01.tif" % (data_dir, plate, file_name), plugin="tifffile")
        img_nuclear_seg = skio.imread("%s%s/seg_tif/%s_seg.tif" % (data_dir1, plate, file_name), plugin="tifffile").astype(int)

        nuclear_props = regionprops(img_nuclear_seg, img_hoechst)
        MYC_props = regionprops(img_nuclear_seg, img_MYC)

        centroid_nuclear = [nuclear_props[x].centroid for x in range(len(nuclear_props))]
        area_nuclear = [nuclear_props[x].area for x in range(len(nuclear_props))]
        perimeter_nuclear = [nuclear_props[x].perimeter for x in range(len(nuclear_props))]
        mean_int_nuclear = [nuclear_props[x].intensity_mean for x in range(len(nuclear_props))]
        mean_int_IF = [MYC_props[x].intensity_mean for x in range(len(MYC_props))]
        circ_nuclear = [(4 * math.pi * area_nuclear[x]) / (perimeter_nuclear[x] ** 2) for x in range(len(area_nuclear))]

        data_temp = pd.DataFrame({'plate': [plate] * len(area_nuclear),
                                  'well': ['%s%s' % (row, column)] * len(area_nuclear),
                                  'nuclear': range(len(area_nuclear)), 'FOV': [j+1] * len(area_nuclear),
                                  'n_nuclear_convex_dilation': [0] * len(area_nuclear),
                                  'centroid_nuclear': centroid_nuclear, 'area_nuclear': area_nuclear,
                                  'circ_nuclear': circ_nuclear, 'mean_int_nuclear': mean_int_nuclear,
                                  'mean_int_MYC': mean_int_IF})
        data = pd.concat([data, data_temp], axis=0)

data.to_csv('%s%s.txt' % (output_dir, plate), index=False, sep='\t')
print("DONE!")
