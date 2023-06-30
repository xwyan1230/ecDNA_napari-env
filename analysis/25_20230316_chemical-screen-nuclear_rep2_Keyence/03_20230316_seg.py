import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import shared.objects as obj
import tifffile as tif
from skimage.morphology import erosion, disk
from skimage.measure import label, regionprops
import math
import shared.segmentation as seg
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230316_analysis_chemical-screen-nuclear_rep2_Keyence/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

# set parameters
otsu_factor = 1.2
circ_threshold = 0.8
min_size = 400
max_size = 1800

plate = 'HSR_2hr'
num_total_samples = 240
start_num = 0
for f in range(num_total_samples):
    fov = start_num + f
    file_name = 'Image_XY0%s' % (fov+1) if fov<9 else 'Image_XY%s' % (fov+1)
    print(file_name)
    img_hoechst = skio.imread("%s%s/raw_img/%s_CH1.tif" % (data_dir, plate, file_name), plugin="tifffile")[:, :, 2] # 1440x1920
    img_MYC = skio.imread("%s%s/raw_img/%s_CH3.tif" % (data_dir, plate, file_name), plugin="tifffile")[:, :, 0]

    # nuclear seg
    img_nuclear_seg = seg.cell_seg_fluorescent(img_hoechst, otsu_factor=otsu_factor, max_size=max_size, min_size=min_size, circ_thresh=circ_threshold).astype(int)

    if not os.path.exists("%s%s/seg_tif/" % (output_dir, plate)):
        os.makedirs("%s%s/seg_tif/" % (output_dir, plate))
    tif.imwrite("%s%s/seg_tif/%s_seg.tif" % (output_dir, plate, file_name), img_nuclear_seg)

    if not os.path.exists("%s%s/color_img/" % (output_dir, plate)):
        os.makedirs("%s%s/color_img/" % (output_dir, plate))

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 255])
    viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 255])
    viewer.add_image(img_nuclear_seg, blending='additive')
    plt.imsave("%s%s/color_img/%s_seg.tiff" % (output_dir, plate, file_name), dis.blending(viewer))
    viewer.close()
    # napari.run()
