from skimage.measure import label, regionprops
from skimage.morphology import medial_axis, dilation, erosion, binary_erosion, binary_dilation, extrema, disk
import pandas as pd
import numpy as np
import shared.image as ima
import skimage.io as skio
import shared.dataframe as dat
import random
import shared.segmentation as seg
import tifffile as tif
import shared.math as mat
from skimage import segmentation
import shared.display as dis
from skimage.filters import threshold_otsu, threshold_local, sobel
import shared.objects as obj
import matplotlib.pyplot as plt
import os
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221214_analysis_Natasha_GBMec_dual-funnel/221017_GBMEC_3hrtreatment_EGFR_interphaseFISH/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder
sample = 'slide1_500nMJQ1'
start_fov = 1
total_fov = 50

# cell info
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, Paul scope 2048x2048:102)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
z_size = 500  # nm (Paul scope)
local_size = 125
rmax = 100
size = 20*local_size + 22

# LOAD Z FILE
data_z = pd.read_csv('%s%s_z.txt' % (data_dir1, sample), na_values=['.'], sep='\t')
data_z['centroid_nuclear'] = [dat.str_to_float(data_z['centroid_nuclear'][i]) for i in range(len(data_z))]

fovs = list(set(data_z['FOV'].tolist()))
total_fov = len(fovs)

# create grid images
img_count = 0
cell_count = 0
img_grid_nuclear = np.zeros(shape=(size, size), dtype=np.uint16)
img_grid_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)
img_grid_nuclear_seg = np.zeros(shape=(size, size), dtype=np.uint16)


def grid_location(num: int):
    row = int(num / 10)
    column = num % 10
    x = row * (2 * local_size + 2) + 2 + local_size
    y = column * (2 * local_size + 2) + 2 + local_size
    return x, y


for f in range(total_fov):
    print("Analyzing %s, analyzing fov %s/%s" % (sample, f+1, total_fov))
    fov = fovs[f]
    data_z_fov = data_z[data_z['FOV'] == fov].copy().reset_index()

    # load images
    file_name = 'TileScan 1_Position %s_RAW' % fov
    im_z_stack_nuclear = skio.imread("%s221017_GBMEC_%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
    im_z_stack_DNAFISH = skio.imread("%s221017_GBMEC_%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
    im_z_stack_seg_convex = skio.imread('%s%s/seg_tif/%s_%s_seg.tif' % (data_dir1, sample, sample, fov),
                                        plugin="tifffile")

    for i in range(len(data_z_fov)):
        print("Analyzing %s, analyzing fov %s, nucleus %s/%s" % (sample, f+1, i+1, len(data_z_fov)))
        z_current = data_z_fov['z'][i]
        label_nuclear = data_z_fov['label_nuclear'][i]
        original_centroid_nuclear = data_z_fov['centroid_nuclear'][i]

        # get images for given z
        img_nuclear_seg_convex = im_z_stack_seg_convex[z_current]
        img_nuclear = im_z_stack_nuclear[z_current]
        img_DNAFISH = im_z_stack_DNAFISH[z_current]

        # get local images
        position = ima.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear, local_size)
        local_nuclear_seg_convex = ima.img_local_seg(img_nuclear_seg_convex, position, label_nuclear)
        local_nuclear = img_nuclear.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH.copy()
        local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
        local_nuclear_props = regionprops(label(local_nuclear_seg_convex))
        local_nuclear_centroid = local_nuclear_props[0].centroid

        if cell_count == 100:
            if not os.path.exists("%s%s/grid/" % (output_dir, sample)):
                os.makedirs("%s%s/grid/" % (output_dir, sample))
            tif.imwrite("%s%s/grid/%s_nuclear_%s.tif" % (output_dir, sample, sample, img_count), img_grid_nuclear)
            tif.imwrite("%s%s/grid/%s_DNAFISH_%s.tif" % (output_dir, sample, sample, img_count), img_grid_DNAFISH)
            tif.imwrite("%s%s/grid/%s_seg_%s.tif" % (output_dir, sample, sample, img_count), img_grid_nuclear_seg)

            viewer = napari.Viewer()
            viewer.add_image(img_grid_nuclear, blending='additive', colormap='blue',
                             contrast_limits=[0, img_grid_nuclear.max()])
            viewer.add_image(img_grid_DNAFISH, blending='additive', colormap='green',
                             contrast_limits=[0, img_grid_DNAFISH.max()])
            plt.imsave('%s%s/grid/%s_%s.tiff' % (output_dir, sample, sample, img_count), dis.blending(viewer))
            viewer.close()

            cell_count = 0
            img_count += 1
            img_grid_nuclear = np.zeros(shape=(size, size), dtype=np.uint16)
            img_grid_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)
            img_grid_nuclear_seg = np.zeros(shape=(size, size), dtype=np.uint16)

        elif cell_count < 100:
            paste_to_centroid = [grid_location(cell_count)[0], grid_location(cell_count)[1]]
            print(paste_to_centroid)
            cell_count += 1
            img_grid_nuclear = ima.image_paste_to(img_grid_nuclear, local_nuclear,
                                                  [int(paste_to_centroid[0] - local_nuclear_centroid[0]),
                                                   int(paste_to_centroid[1] - local_nuclear_centroid[1])])
            img_grid_DNAFISH = ima.image_paste_to(img_grid_DNAFISH, local_DNAFISH,
                                                  [int(paste_to_centroid[0] - local_nuclear_centroid[0]),
                                                   int(paste_to_centroid[1] - local_nuclear_centroid[1])])
            img_grid_nuclear_seg = ima.image_paste_to(img_grid_nuclear_seg, local_nuclear_seg_convex,
                                                  [int(paste_to_centroid[0] - local_nuclear_centroid[0]),
                                                   int(paste_to_centroid[1] - local_nuclear_centroid[1])])

if not os.path.exists("%s%s/grid/" % (output_dir, sample)):
    os.makedirs("%s%s/grid/" % (output_dir, sample))
tif.imwrite("%s%s/grid/%s_nuclear_%s.tif" % (output_dir, sample, sample, img_count), img_grid_nuclear)
tif.imwrite("%s%s/grid/%s_DNAFISH_%s.tif" % (output_dir, sample, sample, img_count), img_grid_DNAFISH)
tif.imwrite("%s%s/grid/%s_seg_%s.tif" % (output_dir, sample, sample, img_count), img_grid_nuclear_seg)

viewer = napari.Viewer()
viewer.add_image(img_grid_nuclear, blending='additive', colormap='blue',
                 contrast_limits=[0, img_grid_nuclear.max()])
viewer.add_image(img_grid_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_grid_DNAFISH.max()])
plt.imsave('%s%s/grid/%s_%s.tiff' % (output_dir, sample, sample, img_count), dis.blending(viewer))
viewer.close()

print("DONE!")
