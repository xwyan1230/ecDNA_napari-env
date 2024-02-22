import skimage.io as skio
import napari
import imutils
import shared.image as ima
import shared.display as dis
from skimage.morphology import disk, dilation
import tifffile as tif
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import numpy as np
import pandas as pd
# from napari_animation import Animation
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230706_analysis_radial/"
data_dir = "%sdata/raw/" % master_folder
data_dir1 = "%sdata/seg_tif/" % master_folder
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sdata/" % master_folder

# sample = 'Colo320DM_acidFISH_lamin_3d'
sample = 'Colo320HSR_acidFISH_lamin_3d'
# prefix = '20230614_acidFISH_lamin_ColoDM_DM'
prefix = '20230614_acidFISH_lamin_ColoHSR_HSR'
total_fov = 6
start_fov = 1

df = pd.read_csv('%s%s/%s_centroid_chosen.txt' % (data_dir2, sample, sample), na_values=['.'], sep='\t')

n_nuclear_convex_dilation = 2
local_size = 200
size = 20*local_size + 22


def grid_location(num: int):
    row = int(num / 10)
    column = num % 10
    x = row * (2 * local_size + 2) + 2 + local_size
    y = column * (2 * local_size + 2) + 2 + local_size
    return x, y


img_count = 0
count = 0
img_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
img_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)
img_laminB = np.zeros(shape=(size, size), dtype=np.uint16)

for i in range(len(df)):
    if count == 100:
        if not os.path.exists("%sgrid/%s/" % (output_dir, sample)):
            os.makedirs("%sgrid/%s/" % (output_dir, sample))
        tif.imwrite("%s/grid/%s/%s_hoechst_%s.tif" % (output_dir, sample, sample, img_count), img_hoechst)
        tif.imwrite("%s/grid/%s/%s_DNAFISH_%s.tif" % (output_dir, sample, sample, img_count), img_DNAFISH)
        tif.imwrite("%s/grid/%s/%s_laminB_%s.tif" % (output_dir, sample, sample, img_count), img_laminB)

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
        viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
        plt.imsave('%sgrid/%s/%s_%s.tiff' % (output_dir, sample, sample, img_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
        viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
        viewer.add_image(img_laminB, blending='additive', colormap='red', contrast_limits=[0, img_laminB.max()])
        plt.imsave('%sgrid/%s/%s_%s_withLaminB.tiff' % (output_dir, sample, sample, img_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
        plt.imsave('%sgrid/%s/%s_DNAFISH_%s.tiff' % (output_dir, sample, sample, img_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
        plt.imsave('%sgrid/%s/%s_hoechst_%s.tiff' % (output_dir, sample, sample, img_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_laminB, blending='additive', colormap='red', contrast_limits=[0, img_laminB.max()])
        plt.imsave('%sgrid/%s/%s_laminB_%s.tiff' % (output_dir, sample, sample, img_count), dis.blending(viewer))
        viewer.close()

        count = 0
        img_count += 1
        img_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
        img_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)
        img_laminB = np.zeros(shape=(size, size), dtype=np.uint16)

    fov = int(df['FOV'][i])
    z = int(df['z'][i])
    print('fov:%s z:%s' % (fov, z))
    if z < 10:
        file_name = '%s_%s_z0%s' % (prefix, fov, z)
    else:
        file_name = '%s_%s_z%s' % (prefix, fov, z)
    img_nuclear_bgc = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_DNAFISH_bgc = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_laminB_bgc = skio.imread("%s%s/%s_ch02.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_seg = skio.imread("%s%s/%s_seg.tif" % (data_dir1, sample, file_name), plugin="tifffile")

    original_centroid_nuclear = [df['centroid_x'][i], df['centroid_y'][i]]
    position = ima.img_local_position(img_seg, original_centroid_nuclear, local_size)
    local_nuclear_seg_convex = ima.img_local_seg(img_seg, position, df['label'][i])
    local_nuclear = img_nuclear_bgc.copy()[position[0]:position[1], position[2]:position[3]]
    local_DNAFISH = img_DNAFISH_bgc.copy()[position[0]:position[1], position[2]:position[3]]
    local_laminB = img_laminB_bgc.copy()[position[0]:position[1], position[2]:position[3]]
    local_nuclear_props = regionprops(label(local_nuclear_seg_convex))
    local_centroid = local_nuclear_props[0].centroid
    local_seg_filter = dilation(local_nuclear_seg_convex, disk(10))
    local_nuclear[local_seg_filter == 0] = 0
    local_DNAFISH[local_seg_filter == 0] = 0
    local_laminB[local_seg_filter == 0] = 0

    if count < 100:
        paste_to_centroid = [grid_location(count)[0], grid_location(count)[1]]
        print(paste_to_centroid)
        count += 1
        img_hoechst = ima.image_paste_to(img_hoechst, local_nuclear,
                                         [int(paste_to_centroid[0] - local_centroid[0]),
                                          int(paste_to_centroid[1] - local_centroid[1])])
        img_DNAFISH = ima.image_paste_to(img_DNAFISH, local_DNAFISH,
                                         [int(paste_to_centroid[0] - local_centroid[0]),
                                                      int(paste_to_centroid[1] - local_centroid[1])])
        img_laminB = ima.image_paste_to(img_laminB, local_laminB,
                                         [int(paste_to_centroid[0] - local_centroid[0]),
                                          int(paste_to_centroid[1] - local_centroid[1])])

if not os.path.exists("%sgrid/%s/" % (output_dir, sample)):
    os.makedirs("%sgrid/%s/" % (output_dir, sample))
tif.imwrite("%s/grid/%s/%s_hoechst_%s.tif" % (output_dir, sample, sample, img_count), img_hoechst)
tif.imwrite("%s/grid/%s/%s_DNAFISH_%s.tif" % (output_dir, sample, sample, img_count), img_DNAFISH)
tif.imwrite("%s/grid/%s/%s_laminB_%s.tif" % (output_dir, sample, sample, img_count), img_laminB)

viewer = napari.Viewer()
viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
plt.imsave('%sgrid/%s/%s_%s.tiff' % (output_dir, sample, sample, img_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
viewer.add_image(img_laminB, blending='additive', colormap='red', contrast_limits=[0, img_laminB.max()])
plt.imsave('%sgrid/%s/%s_%s_withLaminB.tiff' % (output_dir, sample, sample, img_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
plt.imsave('%sgrid/%s/%s_DNAFISH_%s.tiff' % (output_dir, sample, sample, img_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
plt.imsave('%sgrid/%s/%s_hoechst_%s.tiff' % (output_dir, sample, sample, img_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_laminB, blending='additive', colormap='red', contrast_limits=[0, img_laminB.max()])
plt.imsave('%sgrid/%s/%s_laminB_%s.tiff' % (output_dir, sample, sample, img_count), dis.blending(viewer))
viewer.close()

print("DONE!")