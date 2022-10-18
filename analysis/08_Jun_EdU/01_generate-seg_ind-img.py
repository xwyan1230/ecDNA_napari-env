import shared.segmentation as seg
from skimage.morphology import dilation
import matplotlib.pyplot as plt
import math
import skimage.io as skio
import tifffile as tif
import shared.display as dis
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221017_periphery-localization_analysis/20220927_Jun_GBM39_EGFR/"
sample = 'gbm39ec con'
master_path = '%s%s/' % (master_folder, sample)
end_fov = 10
start_fov = 1
total_fov = end_fov - start_fov + 1
save_path = master_path
# cell info
pixel_size = 160  # nm (sp8 confocal 3144x3144:58.7, Paul scope 2048x2048:102, guess for confocal 40x)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
# z_size = 500  # nm (Paul scope)
# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
n_nuclear_convex_dilation = 3
local_size = 150

# IMAGING ANALYSIS
for f in range(total_fov):
    fov = f + start_fov
    print("Analyzing %s, start nuclear segmentation FOV %s/%s" % (sample, fov, total_fov))
    file_prefix = "40x %s-%s" % (sample, fov)
    # LOAD IMAGE
    img_nuclear = skio.imread("%sC2-%s.tif" % (master_path, file_prefix), plugin="tifffile")
    img_DNAFISH = skio.imread("%sC1-%s.tif" % (master_path, file_prefix), plugin="tifffile")
    img_RPA = skio.imread("%sC3-%s.tif" % (master_path, file_prefix), plugin="tifffile")
    img_EdU = skio.imread("%sC4-%s.tif" % (master_path, file_prefix), plugin="tifffile")

    # Nuclear segmentation
    img_nuclear_seg = seg.nuclear_seg(img_nuclear, local_factor=799,
                                      min_size=min_size_nuclear, max_size=3*max_size_nuclear)
    img_nuclear_seg_convex = seg.obj_to_convex(img_nuclear_seg)
    # img_nuclear_seg_convex = dilation(img_nuclear_seg_convex)

    # DNAFISH segmentation
    img_DNAFISH_seg = seg.puncta_seg(img_DNAFISH, img_nuclear_seg_convex, local_size)

    # RPA segmentation
    """img_RPA_seg = seg.puncta_seg(img_RPA, img_nuclear_seg_convex, local_size, local_cycle_factor=7,
                                 min_threshold_factor=2.75, min_size=12, min_threshold_first_round_factor=1)"""

    viewer = napari.view_image(img_nuclear, blending='additive', colormap='blue')
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green')
    # viewer.add_image(img_IF, blending='additive', colormap='magenta')
    plt.imsave('%s/%s_%s_napari-img.tiff' % (save_path, file_prefix, fov), dis.blending(viewer))
    viewer.close()

    viewer1 = napari.view_image(img_nuclear_seg_convex, blending='additive', colormap='blue')
    viewer1.add_image(img_DNAFISH_seg, blending='additive', colormap='green', contrast_limits=[0, 1])
    plt.imsave('%s/%s_%s_napari-seg.tiff' % (save_path, file_prefix, fov), dis.blending(viewer1))
    viewer1.close()

    """viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue')
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green')
    viewer.add_image(img_RPA, blending='additive', colormap='red')
    viewer.add_image(img_EdU, blending='additive', colormap='magenta')
    viewer.add_image(img_RPA_seg, blending='additive', contrast_limits=[0, 1])
    napari.run()"""

    tif.imwrite('%s%s_nuclear_seg.tif' % (master_path, file_prefix), img_nuclear_seg_convex)
    tif.imwrite('%s%s_DNAFISH_seg.tif' % (master_path, file_prefix), img_DNAFISH_seg)
    # tif.imwrite('%s%s_RPA_seg.tif' % (master_path, file_prefix), img_RPA_seg)

print("DONE!")
