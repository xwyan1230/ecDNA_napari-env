import shared.segmentation as seg
from skimage.morphology import disk, dilation
import math
import numpy as np
import skimage.io as skio
import tifffile as tif
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220325_Natasha_THZ1/"
sample = 'DMSO'
subfolder = '01_raw'
savefolder = '02_seg'
total_fov = 20
total_z = 22
# cell info
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, Paul scope 2048x2048:102)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
z_size = 500 # nm (Paul scope)
# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
n_nuclear_convex_dilation = 3

# IMAGING ANALYSIS
for fov in range(total_fov):
    print("Start nuclear segmentation FOV %s/%s" % (fov + 1, total_fov))
    # LOAD IMAGE
    img_nuclear = skio.imread("%s%s/%s/%s_RAW_ch00_fov%s.tif" %
                              (master_folder, sample, subfolder, sample, fov), plugin="tifffile")
    img_DNAFISH = skio.imread("%s%s/%s/%s_RAW_ch01_fov%s.tif" %
                              (master_folder, sample, subfolder, sample, fov), plugin="tifffile")

    # Perform nuclear segmentation
    nuclear_seg_convex = np.zeros(shape=(img_nuclear.shape[0], img_nuclear.shape[1], img_nuclear.shape[2]),
                                  dtype=np.uint16)
    for z in range(total_z):
        img_nuclear_seg = seg.nuclear_seg(img_nuclear[z], local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                          max_size=max_size_nuclear)
        img_nuclear_seg_convex = seg.obj_to_convex(img_nuclear_seg)
        img_nuclear_seg_convex = dilation(img_nuclear_seg_convex, disk(n_nuclear_convex_dilation))
        nuclear_seg_convex[z] = img_nuclear_seg_convex

    tif.imsave('%s%s/%s/%s_seg_fov%s.tif' % (master_folder, sample, savefolder, sample, fov), nuclear_seg_convex)

print("DONE!")