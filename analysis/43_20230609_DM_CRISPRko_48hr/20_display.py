import skimage.io as skio
import pandas as pd
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import tifffile as tif
import shared.dataframe as dat
import matplotlib.pyplot as plt
import shared.display as dis
import napari
import math
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sdata/" % master_folder

row = 'D'
sample = 'D11'
hue_order = ['GFP', 'mCherry']
df = pd.read_csv('%s/%s/%s_n4_simplified.txt' % (data_dir1, sample, sample), na_values=['.'], sep='\t')
df['r'] = np.sqrt(df['area_nuclear']/math.pi)
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA']/math.pi)
df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt']/df['r']
# df_sort = df[df['total_area_ecDNA_sqrt_normalized'] > 0.2].copy().reset_index(drop=True)
df_sort = df
df_sample = df_sort[df_sort['group'].isin(hue_order)].copy().reset_index(drop=True)
df = df_sample

feature = ['centroid_nuclear']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

n_nuclear_convex_dilation = 2
local_size = 200
size = 20*local_size + 22


def grid_location(num: int):
    row = int(num / 10)
    column = num % 10
    x = row * (2 * local_size + 2) + 2 + local_size
    y = column * (2 * local_size + 2) + 2 + local_size
    return x, y


img_mCh_count = 0
img_GFP_count = 0
mCh_count = 0
GFP_count = 0
img_mCherry_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
img_mCherry_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)
img_GFP_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
img_GFP_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)

for i in range(len(df)):
    if mCh_count == 100:
        if not os.path.exists("%sgrid/%s/" % (output_dir, sample)):
            os.makedirs("%sgrid/%s/" % (output_dir, sample))
        tif.imwrite("%s/grid/%s/%s_mCherry_hoechst_%s.tif" % (output_dir, sample, sample, img_mCh_count), img_mCherry_hoechst)
        tif.imwrite("%s/grid/%s/%s_mCherry_DNAFISH_%s.tif" % (output_dir, sample, sample, img_mCh_count), img_mCherry_DNAFISH)

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_hoechst, blending='additive', colormap='blue',
                         contrast_limits=[0, max(img_mCherry_hoechst.max(), img_GFP_hoechst.max())])
        viewer.add_image(img_mCherry_DNAFISH, blending='additive', colormap='green',
                         contrast_limits=[0, max(img_mCherry_DNAFISH.max(), img_GFP_DNAFISH.max())])
        plt.imsave('%sgrid/%s/%s_mCherry_%s.tiff' % (output_dir, sample, sample, img_mCh_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_DNAFISH, blending='additive', colormap='green',
                         contrast_limits=[0, max(img_mCherry_DNAFISH.max(), img_GFP_DNAFISH.max())])
        plt.imsave('%sgrid/%s/%s_mCherry_DNAFISH_%s.tiff' % (output_dir, sample, sample, img_mCh_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_hoechst, blending='additive', colormap='blue',
                         contrast_limits=[0, max(img_mCherry_hoechst.max(), img_GFP_hoechst.max())])
        plt.imsave('%sgrid/%s/%s_mCherry_hoechst_%s.tiff' % (output_dir, sample, sample, img_mCh_count), dis.blending(viewer))
        viewer.close()

        mCh_count = 0
        img_mCh_count += 1
        img_mCherry_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
        img_mCherry_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)

    if GFP_count == 100:
        if not os.path.exists("%sgrid/%s/" % (output_dir, sample)):
            os.makedirs("%sgrid/%s/" % (output_dir, sample))
        tif.imwrite("%sgrid/%s/%s_GFP_hoechst_%s.tif" % (output_dir, sample, sample, img_GFP_count), img_GFP_hoechst)
        tif.imwrite("%sgrid/%s/%s_GFP_DNAFISH_%s.tif" % (output_dir, sample, sample, img_GFP_count), img_GFP_DNAFISH)

        viewer = napari.Viewer()
        viewer.add_image(img_GFP_hoechst, blending='additive', colormap='blue',
                         contrast_limits=[0, max(img_mCherry_hoechst.max(), img_GFP_hoechst.max())])
        viewer.add_image(img_GFP_DNAFISH, blending='additive', colormap='green',
                         contrast_limits=[0, max(img_mCherry_DNAFISH.max(), img_GFP_DNAFISH.max())])
        plt.imsave('%sgrid/%s/%s_GFP_%s.tiff' % (output_dir, sample, sample, img_GFP_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_GFP_DNAFISH, blending='additive', colormap='green',
                         contrast_limits=[0, max(img_mCherry_DNAFISH.max(), img_GFP_DNAFISH.max())])
        plt.imsave('%sgrid/%s/%s_GFP_DNAFISH_%s.tiff' % (output_dir, sample, sample, img_GFP_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_GFP_hoechst, blending='additive', colormap='blue',
                         contrast_limits=[0, max(img_mCherry_hoechst.max(), img_GFP_hoechst.max())])
        plt.imsave('%sgrid/%s/%s_GFP_hoechst_%s.tiff' % (output_dir, sample, sample, img_GFP_count), dis.blending(viewer))
        viewer.close()

        GFP_count = 0
        img_GFP_count += 1
        img_GFP_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
        img_GFP_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)

    fov = df['FOV'][i]
    print(fov)
    if fov < 10:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_s0%s' % (row, sample, fov)
        # filename = '20230602_DM_CRISPRko_48hr_DNAFISH_%s_%s_s0%s' % (row, sample, fov)
    else:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_s%s' % (row, sample, fov)
        # filename = '20230602_DM_CRISPRko_48hr_DNAFISH_%s_%s_s%s' % (row, sample, fov)
    img_nuclear_bgc = skio.imread("%s/FISH/%s/%s_ch01.tif" % (data_dir, sample, filename), plugin="tifffile")
    img_DNAFISH_bgc = skio.imread("%s/FISH/%s/%s_ch00.tif" % (data_dir, sample, filename), plugin="tifffile")
    img_seg = skio.imread("%s/seg/%s/seg_tif/%s_%s_seg.tif" % (data_dir, sample, sample, fov), plugin="tifffile")

    original_centroid_nuclear = df['centroid_nuclear'][i]
    position = ima.img_local_position(img_seg, original_centroid_nuclear, local_size)
    local_nuclear_seg_convex = ima.img_local_seg(img_seg, position, df['label'][i])
    local_nuclear = img_nuclear_bgc.copy()[position[0]:position[1], position[2]:position[3]]
    local_DNAFISH = img_DNAFISH_bgc.copy()[position[0]:position[1], position[2]:position[3]]
    local_nuclear_props = regionprops(label(local_nuclear_seg_convex))
    local_centroid = local_nuclear_props[0].centroid
    local_seg_filter = dilation(local_nuclear_seg_convex, disk(10))
    local_nuclear[local_seg_filter == 0] = 0
    local_DNAFISH[local_seg_filter == 0] = 0

    if df['group'][i] == 'mCherry':
        if mCh_count < 100:
            paste_to_centroid = [grid_location(mCh_count)[0], grid_location(mCh_count)[1]]
            print(paste_to_centroid)
            mCh_count += 1
            img_mCherry_hoechst = ima.image_paste_to(img_mCherry_hoechst, local_nuclear,
                                                     [int(paste_to_centroid[0] - local_centroid[0]),
                                              int(paste_to_centroid[1] - local_centroid[1])])
            img_mCherry_DNAFISH = ima.image_paste_to(img_mCherry_DNAFISH, local_DNAFISH,
                                                     [int(paste_to_centroid[0] - local_centroid[0]),
                                                          int(paste_to_centroid[1] - local_centroid[1])])
    elif df['group'][i] == 'GFP':
        if GFP_count < 100:
            print('GFP')
            paste_to_centroid = [grid_location(GFP_count)[0], grid_location(GFP_count)[1]]
            print(paste_to_centroid)
            GFP_count += 1
            img_GFP_hoechst = ima.image_paste_to(img_GFP_hoechst, local_nuclear,
                                                 [int(paste_to_centroid[0] - local_centroid[0]),
                                                          int(paste_to_centroid[1] - local_centroid[1])])
            img_GFP_DNAFISH = ima.image_paste_to(img_GFP_DNAFISH, local_DNAFISH,
                                                 [int(paste_to_centroid[0] - local_centroid[0]),
                                                          int(paste_to_centroid[1] - local_centroid[1])])

    """viewer = napari.Viewer()
    viewer.add_image(img_mCherry_pos_hoechst, blending='additive', colormap='blue', contrast_limits=[0, max(img_mCherry_pos_hoechst.max(), img_mCherry_neg_hoechst.max())])
    viewer.add_image(img_mCherry_pos_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, max(img_mCherry_pos_DNAFISH.max(), img_mCherry_neg_DNAFISH.max())])
    # plt.imsave('%s%s/color_img/%s_%s_img_label.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    napari.run()

    viewer1 = napari.Viewer()
    viewer1.add_image(img_mCherry_neg_hoechst, blending='additive', colormap='blue', contrast_limits=[0, max(img_mCherry_pos_hoechst.max(), img_mCherry_neg_hoechst.max())])
    viewer1.add_image(img_mCherry_neg_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, max(img_mCherry_pos_DNAFISH.max(), img_mCherry_neg_DNAFISH.max())])
    # plt.imsave('%s%s/color_img/%s_%s_img_label.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    napari.run()"""

if not os.path.exists("%sgrid/%s/" % (output_dir, sample)):
    os.makedirs("%sgrid/%s/" % (output_dir, sample))
tif.imwrite("%sgrid/%s/%s_mCherry_hoechst_%s.tif" % (output_dir, sample, sample, img_mCh_count), img_mCherry_hoechst)
tif.imwrite("%sgrid/%s/%s_mCherry_DNAFISH_%s.tif" % (output_dir, sample, sample, img_mCh_count), img_mCherry_DNAFISH)
tif.imwrite("%sgrid/%s/%s_GFP_hoechst_%s.tif" % (output_dir, sample, sample, img_GFP_count), img_GFP_hoechst)
tif.imwrite("%sgrid/%s/%s_GFP_DNAFISH_%s.tif" % (output_dir, sample, sample, img_GFP_count), img_GFP_DNAFISH)

viewer = napari.Viewer()
viewer.add_image(img_mCherry_hoechst, blending='additive', colormap='blue',
                 contrast_limits=[0, max(img_mCherry_hoechst.max(), img_GFP_hoechst.max())])
viewer.add_image(img_mCherry_DNAFISH, blending='additive', colormap='green',
                 contrast_limits=[0, max(img_mCherry_DNAFISH.max(), img_GFP_DNAFISH.max())])
plt.imsave('%sgrid/%s/%s_mCherry_%s.tiff' % (output_dir, sample, sample, img_mCh_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_DNAFISH, blending='additive', colormap='green',
                 contrast_limits=[0, max(img_mCherry_DNAFISH.max(), img_GFP_DNAFISH.max())])
plt.imsave('%sgrid/%s/%s_mCherry_DNAFISH_%s.tiff' % (output_dir, sample, sample, img_mCh_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_hoechst, blending='additive', colormap='blue',
                 contrast_limits=[0, max(img_mCherry_hoechst.max(), img_GFP_hoechst.max())])
plt.imsave('%sgrid/%s/%s_mCherry_hoechst_%s.tiff' % (output_dir, sample, sample, img_mCh_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_GFP_hoechst, blending='additive', colormap='blue',
                 contrast_limits=[0, max(img_mCherry_hoechst.max(), img_GFP_hoechst.max())])
viewer.add_image(img_GFP_DNAFISH, blending='additive', colormap='green',
                 contrast_limits=[0, max(img_mCherry_DNAFISH.max(), img_GFP_DNAFISH.max())])
plt.imsave('%sgrid/%s/%s_GFP_%s.tiff' % (output_dir, sample, sample, img_GFP_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_GFP_DNAFISH, blending='additive', colormap='green',
                 contrast_limits=[0, max(img_mCherry_DNAFISH.max(), img_GFP_DNAFISH.max())])
plt.imsave('%sgrid/%s/%s_GFP_DNAFISH_%s.tiff' % (output_dir, sample, sample, img_GFP_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_GFP_hoechst, blending='additive', colormap='blue',
                 contrast_limits=[0, max(img_mCherry_hoechst.max(), img_GFP_hoechst.max())])
plt.imsave('%sgrid/%s/%s_GFP_hoechst_%s.tiff' % (output_dir, sample, sample, img_GFP_count), dis.blending(viewer))
viewer.close()