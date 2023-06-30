import skimage.io as skio
import numpy as np
from skimage.measure import regionprops
import math
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230316_analysis_chemical-screen-nuclear_rep2_Keyence/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

data = pd.DataFrame(columns=['exp', 'plate', 'well', 'nuclear', 'FOV', 'n_nuclear_convex_dilation',
                             'centroid_nuclear', 'area_nuclear', 'circ_nuclear', 'mean_int_nuclear', 'mean_int_MYC'])

exp = 'HSR_2hr'
plate = 'HSR_2hr'
num_total_samples = 240
start_num = 0
rows = ['B', 'C', 'D', 'E', 'F', 'G']
columns = ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
for f in range(num_total_samples):
    fov = start_num + f
    row_num = int(fov/(4*len(columns)))
    column_num = int(fov/4 - row_num * len(columns))
    row = rows[row_num]
    column = columns[column_num]
    FOV = fov - (row_num*4*len(columns) + 4*column_num)
    print('%s%s, fov%s' % (row, column, FOV))
    file_name = 'Image_XY0%s' % (fov+1) if fov<9 else 'Image_XY%s' % (fov+1)
    print(file_name)
    img_hoechst = skio.imread("%s%s/raw_img/%s_CH1.tif" % (data_dir, plate, file_name), plugin="tifffile")[:, :, 2] # 1440x1920
    img_MYC = skio.imread("%s%s/raw_img/%s_CH3.tif" % (data_dir, plate, file_name), plugin="tifffile")[:, :, 0]
    img_nuclear_seg = skio.imread("%s%s/seg_tif/%s_seg.tif" % (data_dir, plate, file_name), plugin="tifffile").astype(int)

    nuclear_props = regionprops(img_nuclear_seg, img_hoechst)
    MYC_props = regionprops(img_nuclear_seg, img_MYC)

    centroid_nuclear = [nuclear_props[x].centroid for x in range(len(nuclear_props))]
    area_nuclear = [nuclear_props[x].area for x in range(len(nuclear_props))]
    perimeter_nuclear = [nuclear_props[x].perimeter for x in range(len(nuclear_props))]
    mean_int_nuclear = [nuclear_props[x].intensity_mean for x in range(len(nuclear_props))]
    mean_int_IF = [MYC_props[x].intensity_mean for x in range(len(MYC_props))]
    circ_nuclear = [(4 * math.pi * area_nuclear[x]) / (perimeter_nuclear[x] ** 2) for x in range(len(area_nuclear))]


    data_temp = pd.DataFrame({'exp': [exp] * len(area_nuclear),
                              'plate': [plate] * len(area_nuclear),
                              'well': ['%s%s' % (row, column)] * len(area_nuclear),
                              'nuclear': range(len(area_nuclear)), 'FOV': [FOV+1] * len(area_nuclear),
                              'n_nuclear_convex_dilation': [0] * len(area_nuclear),
                              'centroid_nuclear': centroid_nuclear, 'area_nuclear': area_nuclear,
                              'circ_nuclear': circ_nuclear, 'mean_int_nuclear': mean_int_nuclear,
                              'mean_int_MYC': mean_int_IF})
    data = pd.concat([data, data_temp], axis=0)

data.to_csv('%s%s.txt' % (output_dir, plate), index=False, sep='\t')
print("DONE!")
