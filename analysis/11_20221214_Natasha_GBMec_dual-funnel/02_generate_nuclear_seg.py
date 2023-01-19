import shared.segmentation as seg
from skimage.measure import regionprops
import math
import numpy as np
import skimage.io as skio
import tifffile as tif
import os
import pandas as pd
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221214_analysis_Natasha_GBMec_dual-funnel/221017_GBMEC_3hrtreatment_EGFR_interphaseFISH/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder
sample = 'slide2_1uMtriptolide'
start_fov = 1
total_fov = 50

# cell info
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, Paul scope 2048x2048:102)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
z_size = 500  # nm (Paul scope)

# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
n_nuclear_convex_dilation = 0

data = pd.DataFrame(columns=['FOV',
                             'z',
                             'label_nuclear',
                             'centroid_nuclear',
                             'mean_intensity_MYC_DNAFISH_in_nucleus'])

for f in range(total_fov):
    fov = f + start_fov
    print("Analyzing %s, start nuclear segmentation FOV %s/%s" % (sample, f + 1, total_fov))
    file_name = 'TileScan 1_Position %s_RAW' % fov
    im_z_stack_nuclear = skio.imread("%s221017_GBMEC_%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
    im_z_stack_DNAFISH = skio.imread("%s221017_GBMEC_%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")

    total_z = im_z_stack_nuclear.shape[0]

    # Perform nuclear segmentation
    im_z_stack_nuclear_seg_convex = np.zeros(shape=im_z_stack_nuclear.shape, dtype=np.uint16)
    for z in range(total_z):
        img_nuclear_seg = seg.nuclear_seg(im_z_stack_nuclear[z], local_factor=local_factor_nuclear,
                                          min_size=min_size_nuclear, max_size=max_size_nuclear)
        img_nuclear_seg_convex = seg.obj_to_convex(img_nuclear_seg)
        # img_nuclear_seg_convex = dilation(img_nuclear_seg_convex, disk(n_nuclear_convex_dilation))
        im_z_stack_nuclear_seg_convex[z] = img_nuclear_seg_convex

        nuclear_props = regionprops(img_nuclear_seg_convex, im_z_stack_nuclear[z])
        DNAFISH_props = regionprops(img_nuclear_seg_convex, im_z_stack_DNAFISH[z])

        for i in range(len(nuclear_props)):
            label_nuclear = nuclear_props[i].label
            centroid_nuclear = nuclear_props[i].centroid
            FISH_mean_intensity_nuclear = DNAFISH_props[i].mean_intensity
            data.loc[len(data.index)] = [fov, z, label_nuclear, centroid_nuclear, FISH_mean_intensity_nuclear]

    """viewer = napari.Viewer()
    viewer.add_image(im_z_stack_nuclear, blending='additive', colormap='blue', contrast_limits=[0, im_z_stack_nuclear.max()])
    viewer.add_image(im_z_stack_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, im_z_stack_DNAFISH.max()])
    viewer.add_image(im_z_stack_nuclear_seg_convex, blending='additive', colormap='gray', contrast_limits=[0, 1])
    napari.run()"""
    if not os.path.exists("%s%s/seg_tif/" % (output_dir, sample)):
        os.makedirs("%s%s/seg_tif/" % (output_dir, sample))
    tif.imwrite('%s%s/seg_tif/%s_%s_seg.tif' % (output_dir, sample, sample, fov), im_z_stack_nuclear_seg_convex)

data.to_csv('%s%s_centroids.txt' % (output_dir, sample), index=False, sep='\t')

print("DONE!")
