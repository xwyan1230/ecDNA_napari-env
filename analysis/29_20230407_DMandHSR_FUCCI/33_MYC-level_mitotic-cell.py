import napari
import shared.display as dis
import tifffile as tif
import matplotlib.pyplot as plt
import numpy as np
import cv2
import math
import shared.segmentation as seg
import shared.objects as obj
import pandas as pd
import shared.image as ima
from skimage.measure import label, regionprops_table, regionprops
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk
import os
import skimage.io as skio

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'DM_324pos_merge'
local_size = 1500

img_hoechst = skio.imread("%s%s/DM_hoechst_IF_22000-33601.tif" % (data_dir, sample), plugin="tifffile")
img_nuclear_seg = skio.imread("%s%s/DM_nuclear_red_green_IF_22000-33601_seg_convex.tif" % (data_dir, sample), plugin="tifffile")

df = pd.read_csv('%s%s_cellcycle.txt' % (output_dir, sample), na_values=['.'], sep='\t')

df_M = df[df['cellcycle'] == 'M'].copy().reset_index(drop=True)

mitosis = ['prophase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'prophase', 'prophase', 'NA', 'telophase', 'NA', 'NA', 'NA', 'NA', 'NA', 'prophase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'telophase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'prophase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'prophase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'prophase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'NA', 'NA', 'NA', 'anaphase', 'NA', 'metaphase', 'NA', 'metaphase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'metaphase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'NA', 'NA', 'metaphase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'NA', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'NA', 'NA', 'metaphase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'NA', 'NA', 'NA', 'NA', 'anaphase', 'NA', 'NA', 'NA', 'NA', 'anaphase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'metaphase', 'metaphase', 'NA', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'NA', 'metaphase', 'NA', 'metaphase', 'NA', 'metaphase', 'NA', 'NA', 'metaphase', 'NA', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'anaphase', 'metaphase', 'metaphase', 'NA', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'NA', 'NA', 'NA', 'metaphase', 'metaphase', 'NA', 'metaphase', 'metaphase', 'NA', 'NA', 'metaphase', 'metaphase', 'metaphase', 'NA', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'NA', 'metaphase', 'NA', 'metaphase', 'metaphase', 'NA', 'NA', 'anaphase', 'metaphase', 'NA', 'NA', 'NA', 'metaphase', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'metaphase', 'prophase', 'metaphase', 'metaphase', 'metaphase', 'NA', 'metaphase', 'metaphase', 'metaphase', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'NA', 'prophase', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'metaphase', 'NA', 'metaphase', 'NA', 'metaphase', 'NA', 'prophase', 'metaphase', 'NA', 'anaphase', 'NA', 'NA', 'metaphase', 'metaphase', 'NA', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'NA', 'NA', 'metaphase', 'metaphase', 'NA', 'NA', 'metaphase', 'metaphase', 'metaphase', 'NA', 'telophase', 'metaphase', 'NA', 'NA', 'metaphase', 'NA', 'NA', 'NA', 'metaphase', 'anaphase', 'NA', 'metaphase', 'NA', 'NA', 'NA', 'NA', 'metaphase', 'metaphase', 'metaphase', 'prophase', 'metaphase', 'metaphase', 'telophase', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'NA', 'metaphase', 'metaphase', 'NA', 'metaphase', 'telophase', 'NA', 'metaphase', 'NA', 'metaphase', 'metaphase', 'prophase', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'telophase', 'metaphase', 'metaphase', 'NA', 'metaphase', 'NA', 'metaphase', 'NA', 'metaphase', 'metaphase', 'metaphase', 'NA', 'metaphase', 'NA', 'NA', 'metaphase', 'metaphase', 'metaphase', 'metaphase', 'NA', 'NA', 'NA', 'NA', 'prophase', 'metaphase', 'metaphase', 'NA', 'NA', 'NA', 'NA']

"""start_i = 395
for temp_i in range(len(df_M)-start_i+1):
    i = temp_i + start_i
    print(i)
    print(mitosis)
    label = df_M['label'][i]
    img_temp = np.zeros_like(img_nuclear_seg)
    img_temp[img_nuclear_seg == label] = 1
    props_temp = regionprops(img_temp)
    centroid = np.array(props_temp[0].centroid)
    position = ima.img_local_position(img_nuclear_seg, centroid, local_size)
    local_nuclear_seg = ima.img_local_seg(img_nuclear_seg, position, label)
    local_nuclear = img_hoechst.copy()
    local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
    local_props = regionprops(local_nuclear_seg.astype(int))
    local_centroid = local_props[0].centroid

    viewer = napari.Viewer()
    viewer.add_image(local_nuclear, blending='additive', colormap='blue', contrast_limits=[0, local_nuclear.max()])
    viewer.add_points(local_centroid, face_color='red')
    shapes = viewer.add_shapes(name='Shapes', ndim=2)
    napari.run()
    shapes_layer = viewer.layers['Shapes']

    if len(shapes.data) == 0:
        mitosis.append('NA')
    elif len(shapes.data) == 1:
        mitosis.append('prophase')
    elif len(shapes.data) == 2:
        mitosis.append('metaphase')
    elif len(shapes.data) == 3:
        mitosis.append('anaphase')
    elif len(shapes.data) == 4:
        mitosis.append('telophase')
    else:
        mitosis.append('NA')"""

df_M['mitosis'] = mitosis
df_M.to_csv('%s%s_mitosis.txt' % (output_dir, sample), index=False, sep='\t')
print("DONE!")


