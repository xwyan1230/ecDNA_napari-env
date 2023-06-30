import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import shared.objects as obj
import pandas as pd
import tifffile as tif
from skimage.morphology import erosion, disk
from skimage.measure import label, regionprops
import math
import shared.segmentation as seg
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230508_analysis_chemical-screen-nuclear_24hr/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

# set parameters
pixel_size = 142  # nm (sp8 confocal 40x 2048*2048)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.8, 1.5]  # used to filter nucleus, DM [0.8, 1.3], HSR [0.8, 1.5]
n_nuclear_convex_dilation = 12
convex_conversion_threshold = 0.85
local_size = 200
circ_threshold = 0.8

# segmentation
local_factor_nuclear = 77  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi

plate = 'HSR_24hr'
seq = pd.read_csv('%s%s_sequence.txt' % (data_dir, plate), na_values=['.'], sep='\t')

samples = seq['sample'].tolist()
start_sample = 79
# 24,64,54,45,33

for s in range(20):
    print(s+start_sample)
    sample = samples[s+start_sample]
    print(sample)

    file_name = '20230326_HSR_24hr_chemical-screen-nuclear_oldMYC_%s_RAW' % sample
    # file_name = '20230325_DM_chemical-screen-nuclear_24hr_oldMYC_%s_RAW' % sample
    img_hoechst_stack = skio.imread("%s%s/%s_ch01.tif" % (data_dir, plate, file_name), plugin="tifffile")
    img_MYC_stack = skio.imread("%s%s/%s_ch00.tif" % (data_dir, plate, file_name), plugin="tifffile")

    for fov in range(4):
        img_hoechst = img_hoechst_stack[:, :, fov]
        img_MYC = img_MYC_stack[:, :, fov]
        # img_MYC = np.concatenate([np.zeros(shape=[5, 2048]), img_MYC], axis=0)[:2048, :2048]

        # nuclear seg
        img_nuclear_seg = seg.nuclear_seg1(img_hoechst, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                           max_size=max_size_nuclear)

        img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg,
                                                                           threshold=convex_conversion_threshold))

        nuclear_props = regionprops(img_nuclear_seg_convex)
        img_nuclear_seg_final = np.zeros_like(img_nuclear_seg_convex)
        for k in range(len(nuclear_props)):
            circ_nuclear = (4 * math.pi * nuclear_props[k].area) / (nuclear_props[k].perimeter ** 2)
            if circ_nuclear > circ_threshold:
                img_nuclear_seg_final[img_nuclear_seg_convex == nuclear_props[k].label] = nuclear_props[k].label
        img_nuclear_seg_final = obj.label_resort(img_nuclear_seg_final)

        if not os.path.exists("%s%s/seg_tif/" % (output_dir, plate)):
            os.makedirs("%s%s/seg_tif/" % (output_dir, plate))
        tif.imwrite("%s%s/seg_tif/%s_seg_%s.tif" % (output_dir, plate, file_name, fov), img_nuclear_seg_final)

        if not os.path.exists("%s%s/color_img/" % (output_dir, plate)):
            os.makedirs("%s%s/color_img/" % (output_dir, plate))

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
        plt.imsave("%s%s/color_img/%s_img_%s.tiff" % (output_dir, plate, file_name, fov), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
        viewer.add_image(img_nuclear_seg_final, blending='additive')
        plt.imsave("%s%s/color_img/%s_seg_%s.tiff" % (output_dir, plate, file_name, fov), dis.blending(viewer))
        viewer.close()
