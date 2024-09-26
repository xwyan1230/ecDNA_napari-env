import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import tifffile as tif
from skimage.measure import regionprops, label
from skimage.morphology import binary_dilation, disk, dilation, erosion
import numpy as np
import shared.segmentation as seg
import shared.objects as obj
import pandas as pd
import shared.image as ima
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240227_analysis_GFPandmCherry_mixing-ratio_test/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

folder = '384well'
sample = 'XY01'
fov = 15
manual = 'N'

data = pd.DataFrame()
fov_lst = []
GFP_lst = []
mCherry_lst = []
emiRFP670_lst = []

file_name = 'Image_%s_000%s' % (sample, fov)
img_green = skio.imread("%s%s/%s/%s_CH2.tif" % (data_dir, folder, sample, file_name), plugin="tifffile")[:, :, 1]
img_red = skio.imread("%s%s/%s/%s_CH3.tif" % (data_dir, folder, sample, file_name), plugin="tifffile")[:, :, 0]
img_trans = skio.imread("%s%s/%s/%s_CH4.tif" % (data_dir, folder, sample, file_name), plugin="tifffile")

if manual == 'YES':
    if os.path.exists("%s/%s/manual_seg/%s/%s_m-seg.tif" % (output_dir, folder, sample, file_name)):
        img_seg = skio.imread("%s%s/manual_seg/%s/%s_m-seg.tif" % (output_dir, folder, sample, file_name),
                              plugin="tifffile")
    else:
        img_seg = np.zeros_like(img_trans)

    viewer = napari.Viewer()
    viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, 30000])
    viewer.add_image(img_red, blending='additive', colormap='red', contrast_limits=[0, 30000])
    viewer.add_image(img_trans, blending='additive', contrast_limits=[0, 65535])
    viewer.add_image(img_seg, blending='additive', contrast_limits=[0, 1])
    shapes_add = viewer.add_shapes(name='add', ndim=2)
    napari.run()

    img_seg = erosion(img_seg, disk(6))
    img_seg = ima.napari_add_or_remove_obj(shapes_add.data, 'add', img_seg)
    img_seg = dilation(img_seg, disk(6))
    viewer = napari.Viewer()
    viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, 30000])
    viewer.add_image(img_red, blending='additive', colormap='red', contrast_limits=[0, 30000])
    viewer.add_image(img_trans, blending='additive', contrast_limits=[0, 65535])
    viewer.add_image(img_seg, blending='additive', contrast_limits=[0, 1])
    napari.run()
    if not os.path.exists("%s%s/manual_seg/%s/" % (output_dir, folder, sample)):
        os.makedirs("%s%s/manual_seg/%s/" % (output_dir, folder, sample))
    tif.imwrite("%s/%s/manual_seg/%s/%s_m-seg.tif" % (output_dir, folder, sample, file_name), img_seg)
else:
    img_seg = skio.imread("%s%s/manual_seg/%s/%s_m-seg.tif" % (output_dir, folder, sample, file_name), plugin="tifffile")
label_img = img_seg
green_props = regionprops(label_img, img_green)
red_props = regionprops(label_img, img_red)
fov_lst = fov_lst + [fov] * len(green_props)
GFP_lst = GFP_lst + [green_props[j].intensity_mean for j in range(len(green_props))]
mCherry_lst = mCherry_lst + [red_props[j].intensity_mean for j in range(len(red_props))]

data['fov'] = fov_lst
data['GFP'] = GFP_lst
data['mCherry'] = mCherry_lst

data.to_csv('%s/%s/txt/%s_%s_manual.txt' % (output_dir, folder, sample, fov), index=False, sep='\t')
print("DONE!")
