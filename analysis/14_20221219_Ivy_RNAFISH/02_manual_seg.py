import skimage.io as skio
import napari
import tifffile as tif
import shared.display as dis
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
import matplotlib.pyplot as plt
import shared.image as ima
import numpy as np
import shared.segmentation as seg
import shared.objects as obj
import os
import math

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221219_analysis_Ivy_RNAFISH/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder


def fov_to_str(fov):
    if fov < 10:
        out = '00%s' % fov
    else:
        out = '0%s' % fov
    return out


sample = 'PC3'
# fov_lst = [fov_to_str(i) for i in np.arange(2, 37, 1)] + ['Mn_'+fov_to_str(i) for i in np.arange(1, 7, 1)]  # Colo320DM
# fov_lst = [fov_to_str(i) for i in np.arange(1, 31, 1)]  # Colo320HSR
# fov_lst = ['MycI2_'+fov_to_str(i) for i in np.arange(1, 17, 1)]  # HCT116
fov_lst = list(np.arange(1, 9, 1)) + list(np.arange(10, 13, 1)) + [fov_to_str(i) for i in np.arange(13, 16, 1)] + [fov_to_str(i) for i in np.arange(17, 26, 1)] + [fov_to_str(i) for i in np.arange(27, 36, 1)] # PC3
# fov_lst = [fov_to_str(i) for i in np.arange(1, 2, 1)] + ['Mn_'+fov_to_str(i) for i in np.arange(1, 8, 1)] + ['Mn_'+fov_to_str(i) for i in np.arange(9, 19, 1)]  # PC9

# set parameters
pixel_size = 58.7  # nm (sp8 confocal 3144x3144)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.4, 0.8]  # used to filter nucleus
convex_conversion_threshold = 0.9

# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi

for fov in range(len(fov_lst)):
    print(fov)

    file_name = '%s_%s_Lng_SVCC_Processed001_RAW' % (sample, fov_lst[fov])
    img_hoechst = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_RNAFISH = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_RNAFISH_seg = np.zeros_like(img_RNAFISH)

    img_nuclear_seg = seg.nuclear_seg1(img_hoechst, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                       max_size=max_size_nuclear)
    img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg,
                                                                       threshold=convex_conversion_threshold))

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_RNAFISH, blending='additive', colormap='red', contrast_limits=[0, img_RNAFISH.max()])
    viewer.add_image(img_nuclear_seg_convex, blending='additive')
    shapes_add = viewer.add_shapes(name='nuclear seg add', ndim=2)
    shapes_seg = viewer.add_shapes(name='nuclear seg remove', ndim=2)
    napari.run()

    img_nuclear_seg_convex = ima.napari_add_or_remove_obj(shapes_seg.data, 'remove', img_nuclear_seg_convex)
    img_nuclear_seg_convex = ima.napari_add_or_remove_obj(shapes_add.data, 'add', img_nuclear_seg_convex)

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_RNAFISH, blending='additive', colormap='red', contrast_limits=[0, img_RNAFISH.max()])
    viewer.add_image(img_nuclear_seg_convex, blending='additive')
    points = viewer.add_points(name='Points', ndim=2)
    napari.run()

    for j in points.data:
        img_RNAFISH_seg[int(j[0]), int(j[1])] = 1
    img_RNAFISH_seg = dilation(img_RNAFISH_seg, disk(3))

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_RNAFISH, blending='additive', colormap='red', contrast_limits=[0, img_RNAFISH.max()])
    viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_RNAFISH_seg, blending='additive', colormap='red', contrast_limits=[0, 1])
    napari.run()

    if not os.path.exists("%s%s/seg_tif_manual/" % (output_dir, sample)):
        os.makedirs("%s%s/seg_tif_manual/" % (output_dir, sample))
    tif.imwrite("%s%s/seg_tif_manual/%s_%s_seg.tif" % (output_dir, sample, sample, fov_lst[fov]), img_nuclear_seg_convex)
    tif.imwrite("%s%s/seg_tif_manual/%s_%s_RNAseg.tif" % (output_dir, sample, sample, fov_lst[fov]), img_RNAFISH_seg)

    if not os.path.exists("%s%s/color_img_manual/" % (output_dir, sample)):
        os.makedirs("%s%s/color_img_manual/" % (output_dir, sample))
    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_RNAFISH, blending='additive', colormap='red', contrast_limits=[0, img_RNAFISH.max()])
    plt.imsave('%s%s/color_img_manual/%s_%s_img.tiff' % (output_dir, sample, sample, fov_lst[fov]), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_RNAFISH_seg, blending='additive', colormap='red', contrast_limits=[0, 1])
    plt.imsave('%s%s/color_img_manual/%s_%s_seg.tiff' % (output_dir, sample, sample, fov_lst[fov]), dis.blending(viewer))
    viewer.close()

