import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.dataframe as dat
import shared.math as mat
import napari
import math
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240205_analysis_sp8_dWee1andiWee1/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sprocessed/" % master_folder
output_dir = "%sprocessed/" % master_folder

n_nuclear_convex_dilation = 0
local_size = 200
start_i = 301

samples = ['DMSO_24hr_1', 'dWee1_1uM_24hr_1', 'dWee1_1uM_24hr_2']

df = pd.DataFrame()
for sample in samples:
    df1 = pd.read_csv('%s/txt/%s_n0.txt' % (data_dir1, sample), na_values=['.'], sep='\t')
    df = pd.concat([df, df1], axis=0)

# df_random = df.sample(frac=1)
df_random = pd.read_csv('%s/txt/random.txt' % data_dir1, na_values=['.'], sep='\t')
# df_random.to_csv('%s/txt/random.txt' % output_dir, index=False, sep='\t')
catagory_lst1 = [2, 3, 1, 5, 1, 3, 1, 4, 3, 4, 0, 1, 5, 4, 1, 2, 2, 5, 4, 1, 1, 5, 5, 0, 4, 4, 1, 1, 1, 2, 2, 6, 2, 1,
                 4, 4, 2, 2, 1, 0, 4, 1, 3, 1, 0, 4, 0, 4, 2, 4, 5, 4, 0, 4, 1, 2, 2, 6, 2, 1, 1, 2, 2, 5, 5, 3, 4, 1,
                 6, 0, 1, 1, 0, 1, 1, 0, 6, 0, 2, 6, 1, 5, 1, 1, 6, 3, 4, 0, 4, 0, 2, 3, 6, 1, 3, 4, 1, 4, 2, 6, 5, 2,
                 2, 2, 2, 3, 6, 3, 1, 3, 2, 4, 6, 1, 5, 5, 1, 5, 0, 5, 5, 1, 1, 1, 1, 6, 5, 4, 0, 1, 1, 1, 1, 2, 3, 4,
                 1, 5, 0, 1, 0, 0, 4, 0, 0, 4, 0, 1, 6, 2, 2, 2, 0, 1, 0, 0, 3, 6, 1, 1, 1, 0, 6, 0, 0, 1, 4, 6, 1, 1,
                 6, 1, 1, 0, 0, 2, 1, 1, 1, 0, 0, 0, 0, 1, 3, 1, 0, 0, 0, 6, 6, 5, 5, 2, 0, 1, 1, 4, 3, 0, 3, 1, 1, 1,
                 0, 5, 3, 4, 0, 4, 0, 1, 0, 0, 0, 3, 0, 1, 6, 5, 2, 3, 5, 3, 1, 1, 1, 5, 5, 0, 6, 0, 1, 1, 1, 0, 1, 3,
                 1, 1, 2, 2, 3, 0, 5, 1, 1, 1, 5, 2, 0, 0, 0, 0, 5, 1, 1, 4, 1, 1, 1, 1, 0, 6, 1, 1, 3, 1, 5, 6, 5, 5,
                 4, 1, 1, 1, 1, 1, 1, 5, 1, 1, 3, 4, 1, 0, 0, 1, 1, 3, 3, 6, 1, 6, 1, 0, 6, 1, 2, 0, 1]
print(len(catagory_lst1))

catagory_lst = []
for ii in range(len(df_random)):
    i = start_i + ii
    print(i)
    sample = df_random['sample'].tolist()[i]
    file_name = df_random['file_name'].tolist()[i]
    fov = df_random['FOV'].tolist()[i]
    img_nuclear = skio.imread("%s%s/%s_s%s_ch00.tif" % (data_dir, sample, file_name, fov), plugin="tifffile")
    img_DNAFISH = skio.imread("%s%s/%s_s%s_ch01.tif" % (data_dir, sample, file_name, fov), plugin="tifffile")
    img_seg = skio.imread("%s/seg_tif/%s/%s_%s_seg.tif" % (data_dir1, sample, file_name, fov), plugin="tifffile")
    n = int(df_random['nuclear'].tolist()[i])
    nuclear_props = regionprops(img_seg)
    original_centroid_nuclear = nuclear_props[n].centroid
    position = ima.img_local_position(img_seg, original_centroid_nuclear, local_size)
    local_nuclear_seg = ima.img_local_seg(img_seg, position, nuclear_props[n].label)
    local_nuclear = img_nuclear.copy()
    local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
    local_nuclear[local_nuclear_seg == 0] = 0
    local_DNAFISH = img_DNAFISH.copy()
    local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
    local_DNAFISH[local_nuclear_seg == 0] = 0

    viewer = napari.Viewer()
    # viewer.add_image(local_nuclear, blending='additive', colormap='blue', contrast_limits=[0, 45000])
    viewer.add_image(local_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 25000])
    shapes = viewer.add_shapes(name='Shapes', ndim=2)
    napari.run()
    shapes_layer = viewer.layers['Shapes']
    poly_data = shapes.data
    catagory_lst.append(len(poly_data))
    print(catagory_lst)

df_random['group'] = catagory_lst1 +  catagory_lst
df1 = df_random[df_random['sample'] == 'DMSO_24hr'].copy().reset_index(drop=True)
df2 = df_random[df_random['sample'] == 'dWee1_1uM_24hr'].copy().reset_index(drop=True)

df1.to_csv('%s/txt/DMSO_24hr_group.txt' % output_dir, index=False, sep='\t')
df2.to_csv('%s/txt/dWee1_1uM_24hr_group.txt' % output_dir, index=False, sep='\t')

print("DONE!")