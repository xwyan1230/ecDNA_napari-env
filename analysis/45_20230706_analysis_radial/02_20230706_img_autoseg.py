import skimage.io as skio
import napari
import imutils
import shared.image as ima
import shared.objects as obj
import shared.segmentation as seg
import tifffile as tif
import numpy as np
# from napari_animation import Animation
import os
import math

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230706_analysis_radial/"
data_dir = "%sdata/raw/" % master_folder
output_dir = "%sdata/" % master_folder

# sample = 'Colo320DM_acidFISH_lamin_3d'
sample = 'Colo320HSR_acidFISH_lamin_3d'
# prefix = '20230614_acidFISH_lamin_ColoDM_DM'
prefix = '20230614_acidFISH_lamin_ColoHSR_HSR'
total_fov = 6
start_fov = 1

pixel_size = 56.8  # nm (sp8 confocal 3144x3144)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
convex_conversion_threshold = 0.85
local_size = 200
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000 / (pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000 / (pixel_size * 2)) ** 2 * math.pi

# The microscope reports the following spacing (in µm)
original_spacing = np.array([0.3, 0.05679, 0.05679])
# We downsampled each slice 4x to make the data smaller
rescaled_spacing = original_spacing * [1, 4, 4]
# Normalize the spacing so that pixels are a distance of 1 apart
spacing = rescaled_spacing / rescaled_spacing[2]

lst = [x.split('%s_' % prefix)[1].split('_ch')[0] for x in os.listdir('%s%s' % (data_dir, sample)) if x[-4:] == '.tif']

for f in range(total_fov):
    fov = start_fov + f
    print(fov)
    lst_temp = [int(x.split('_z')[1]) for x in lst if x.split('_z')[0] == '%s' % fov]
    total_z = max(lst_temp) + 1
    for z in range(total_z):
        print(z)
        if z < 10:
            file_name = '%s_%s_z0%s' % (prefix, fov, z)
        else:
            file_name = '%s_%s_z%s' % (prefix, fov, z)
        img_nuclear = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
        img_DNAFISH = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
        img_laminB = skio.imread("%s%s/%s_ch02.tif" % (data_dir, sample, file_name), plugin="tifffile")

        img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                           max_size=max_size_nuclear)
        img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg,
                                                                           threshold=convex_conversion_threshold))

        if not os.path.exists("%sseg_tif/%s/" % (output_dir, sample)):
            os.makedirs("%sseg_tif/%s/" % (output_dir, sample))
        tif.imwrite("%sseg_tif/%s/%s_seg.tif" % (output_dir, sample, file_name), img_nuclear_seg_convex)

        img_DNAFISH_seg = img_DNAFISH > 12000
        img_DNAFISH_seg = obj.remove_small(img_DNAFISH_seg, 6)
        tif.imwrite("%sseg_tif/%s/%s_ecseg1.tif" % (output_dir, sample, file_name), img_DNAFISH_seg)

print("DONE!")

