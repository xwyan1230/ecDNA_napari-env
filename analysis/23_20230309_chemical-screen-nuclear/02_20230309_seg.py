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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230309_analysis_chemical-screen-nuclear/"
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

plate = 'HSR_2hr'
rows = ['G']
columns = ['8', '9', '10', '11']
# columns = ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
num_total_samples = len(rows) * len(columns)
for i in range(num_total_samples):
    row = rows[int(i/len(columns))]
    column = columns[int(i - int(i/len(columns))*len(columns))]
    print("%s: %s%s" % (plate, row, column))
    imgs = [x for x in os.listdir('%s%s/%s/' % (data_dir, plate, row))]
    if '.DS_Store' in imgs:
        imgs.remove('.DS_Store')
    imgs = [x for x in imgs if '%s_%s%s' % (row, row, column) in x]
    n_imgs = int(len(imgs)/2)
    for j in range(n_imgs):
        file_name = '%s_%s%s_%s_RAW' % (row, row, column, j+1)
        img_hoechst = skio.imread("%s%s/%s/%s_ch01.tif" % (data_dir, plate, row, file_name), plugin="tifffile")
        img_MYC = skio.imread("%s%s/%s/%s_ch00.tif" % (data_dir, plate, row, file_name), plugin="tifffile")
        img_hoechst = np.concatenate([np.zeros(shape=[2048, 5]), img_hoechst], axis=1)[:2048, :2048]
        # tif.imwrite("%s%s/%s/%s_ch01_new.tif" % (data_dir, plate, row, file_name), img_hoechst)

        # nuclear seg
        img_nuclear_seg = seg.nuclear_seg1(img_hoechst, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                           max_size=max_size_nuclear)

        img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg,
                                                                           threshold=convex_conversion_threshold))
        img_nuclear_seg_convex = erosion(img_nuclear_seg_convex, disk(3))

        nuclear_props = regionprops(img_nuclear_seg_convex)
        img_nuclear_seg_final = np.zeros_like(img_nuclear_seg_convex)
        for k in range(len(nuclear_props)):
            circ_nuclear = (4 * math.pi * nuclear_props[k].area) / (nuclear_props[k].perimeter ** 2)
            if circ_nuclear > circ_threshold:
                img_nuclear_seg_final[img_nuclear_seg_convex == nuclear_props[k].label] = nuclear_props[k].label
        img_nuclear_seg_final = obj.label_resort(img_nuclear_seg_final)

        if not os.path.exists("%s%s/%s/seg_tif/" % (output_dir, plate, row)):
            os.makedirs("%s%s/%s/seg_tif/" % (output_dir, plate, row))
        tif.imwrite("%s%s/%s/seg_tif/%s_seg_n-3.tif" % (output_dir, plate, row, file_name), img_nuclear_seg_final)

        if not os.path.exists("%s%s/%s/color_img/" % (output_dir, plate, row)):
            os.makedirs("%s%s/%s/color_img/" % (output_dir, plate, row))

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
        plt.imsave("%s%s/%s/color_img/%s_img.tiff" % (output_dir, plate, row, file_name), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
        viewer.add_image(img_nuclear_seg_final, blending='additive')
        plt.imsave("%s%s/%s/color_img/%s_seg.tiff" % (output_dir, plate, row, file_name), dis.blending(viewer))
        viewer.close()
