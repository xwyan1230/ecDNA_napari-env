import shared.segmentation as seg
from skimage.morphology import disk, dilation
from skimage.measure import regionprops
import math
import numpy as np
import skimage.io as skio
import tifffile as tif
import pandas as pd
import os
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220726_BRDfamily_screen/"
sample_lst = ['B11', 'C11', 'D11', 'E11', 'F11', 'G11']
raw_folder = '01_raw'
save_folder = '02_seg'
total_fov = 15
for sample in sample_lst:
    save_path = "%s%s/%s/%s/" % (master_folder, sample[0], sample[1:], save_folder)
    if not os.path.exists(save_path):
        os.makedirs(save_path)
# cell info
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, Paul scope 2048x2048:102)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
z_size = 500  # nm (Paul scope)
# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
n_nuclear_convex_dilation = 3

# IMAGING ANALYSIS
for sample in sample_lst:
    print("Start analyzing %s..." % sample)
    data = pd.DataFrame(columns=['FOV',
                                 'z',
                                 'label_nuclear',
                                 'centroid_nuclear',
                                 'mean_intensity_MYC_DNAFISH_in_nucleus'])
    for fov in range(total_fov):
        print("Start nuclear segmentation FOV %s/%s" % (fov + 1, total_fov))
        # LOAD IMAGE
        im_z_stack_nuclear = skio.imread("%s%s/%s/%s/R%s_RAW_ch00.tif" %
                                         (master_folder, sample[0], sample[1:], raw_folder, fov+1), plugin="tifffile")
        im_z_stack_DNAFISH = skio.imread("%s%s/%s/%s/R%s_RAW_ch01.tif" %
                                         (master_folder, sample[0], sample[1:], raw_folder, fov+1), plugin="tifffile")

        total_z = im_z_stack_nuclear.shape[0]

        # Perform nuclear segmentation
        im_z_stack_nuclear_seg_convex = np.zeros(shape=(im_z_stack_nuclear.shape[0], im_z_stack_nuclear.shape[1],
                                                        im_z_stack_nuclear.shape[2]), dtype=np.uint16)
        for z in range(total_z):
            img_nuclear_seg = seg.nuclear_seg(im_z_stack_nuclear[z], local_factor=local_factor_nuclear,
                                              min_size=min_size_nuclear, max_size=max_size_nuclear)
            img_nuclear_seg_convex = seg.obj_to_convex(img_nuclear_seg)
            img_nuclear_seg_convex = dilation(img_nuclear_seg_convex, disk(n_nuclear_convex_dilation))
            im_z_stack_nuclear_seg_convex[z] = img_nuclear_seg_convex

            nuclear_props = regionprops(img_nuclear_seg_convex, im_z_stack_nuclear[z])
            DNAFISH_props = regionprops(img_nuclear_seg_convex, im_z_stack_DNAFISH[z])

            for i in range(len(nuclear_props)):
                label_nuclear = nuclear_props[i].label
                centroid_nuclear = nuclear_props[i].centroid
                FISH_mean_intensity_nuclear = DNAFISH_props[i].mean_intensity
                data.loc[len(data.index)] = [fov, z, label_nuclear, centroid_nuclear, FISH_mean_intensity_nuclear]

        tif.imwrite(
            "%s%s/%s/%s/R%s_seg.tif" % (master_folder, sample[0], sample[1:], save_folder, fov+1),
            im_z_stack_nuclear_seg_convex)

    data.to_csv('%s%s/%s/%s_centroids.txt' % (master_folder, sample[0], sample[1:], sample), index=False, sep='\t')

print("DONE!")
