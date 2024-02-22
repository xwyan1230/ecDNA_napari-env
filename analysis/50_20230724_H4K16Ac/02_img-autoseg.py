import skimage.io as skio
import napari
import tifffile as tif
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import os
import math
import shared.segmentation as seg
import shared.objects as obj

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230724_analysis_H4K16Ac/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'HSR'
total_fov = 49
start_fov = 0

# set parameters
pixel_size = 58.7  # nm (sp8 confocal 3144x3144)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 0.9]  # used to filter nucleus
n_nuclear_convex_dilation = 2
convex_conversion_threshold = 0.8
local_size = 200

# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi

for f in range(total_fov):
    fov = start_fov + f
    if fov < 10:
        file_name = '20230721_H4K16Ac_%s_H4K16Ac_s0%s' % (sample, fov)
    else:
        file_name = '20230721_H4K16Ac_%s_H4K16Ac_s%s' % (sample, fov)
    img_nuclear = skio.imread("%sraw/%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_DNAFISH = skio.imread("%sraw/%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_IF = skio.imread("%sraw/%s/%s_ch02.tif" % (data_dir, sample, file_name), plugin="tifffile")

    # nuclear segmentation
    img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                       max_size=max_size_nuclear)
    img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg,
                                                                       threshold=convex_conversion_threshold))

    """viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    viewer.add_image(img_nuclear_seg_convex, blending='additive')
    napari.run()"""

    if not os.path.exists("%sseg_tif_new/%s/" % (output_dir, sample)):
        os.makedirs("%sseg_tif_new/%s/" % (output_dir, sample))
    tif.imwrite("%sseg_tif_new/%s/%s_%s_seg.tif" % (output_dir, sample, sample, fov), img_nuclear_seg_convex)

    # ecDNA segmentation
    img_DNAFISH_seg1 = np.zeros_like(img_DNAFISH)
    img_DNAFISH_seg1[img_DNAFISH > 18000] = 1
    img_DNAFISH_seg1 = obj.remove_small(img_DNAFISH_seg1, 10)
    tif.imwrite("%sseg_tif_new/%s/%s_%s_ecseg1.tif" % (output_dir, sample, sample, fov), img_DNAFISH_seg1)

    """viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    viewer.add_image(img_nuclear_seg_convex, blending='additive')
    viewer.add_image(img_DNAFISH_seg1, blending='additive')
    napari.run()"""

    if not os.path.exists("%scolor_img_new/%s/" % (output_dir, sample)):
        os.makedirs("%scolor_img_new/%s/" % (output_dir, sample))
    viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear.max()])
    # viewer.add_image(img_mCherry, blending='additive', colormap='red', contrast_limits=[0, img_mCherry.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    plt.imsave('%scolor_img_new/%s/%s_%s_img.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear.max()])
    viewer.add_image(img_IF, blending='additive', colormap='red', contrast_limits=[0, img_IF.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    plt.imsave('%scolor_img_new/%s/%s_%s_img_with_IF.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_DNAFISH_seg1, blending='additive', colormap='green', contrast_limits=[0, 1])
    plt.imsave('%scolor_img_new/%s/%s_%s_seg.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

print("DONE!")