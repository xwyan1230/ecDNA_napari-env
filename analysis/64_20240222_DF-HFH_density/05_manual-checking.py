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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240222_analysis_DF-HFH_density/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

# set parameters
otsu_factor = 1.1
circ_threshold = 0.5
min_size = 300
max_size = 2000


def img_crop(img, i):
    if i == 0:
        out = img[0:1440, 0:1920]
    elif (i == 1) | (i == 2):
        out = img[0:1440, 576:1920]
    elif (i == 3) | (i == 6):
        out = img[432:1440, 0:1920]
    elif (i == 4) | (i == 5):
        out = img[432:1440, 0:1344]
    elif (i == 7) | (i == 8):
        out = img[432:1440, 576:1920]
    return out


folder = '48hr_384well'
sample_name = 'DMSO_4x8000-rep8'
sample = 'XY15'
n_img = 9
manual = 'YES'

data = pd.DataFrame()
fov_lst = []
hoechst_lst = []
h23_lst = []
azaleaB5_lst = []
emiRFP670_lst = []

for i in range(n_img):
    print(i)
    file_name = 'Image_%s_0000%s' % (sample, i+1)
    img_hoechst = img_crop(skio.imread("%s%s/%s/%s_CH1.tif" % (data_dir, folder, sample, file_name), plugin="tifffile")[:, :, 2], i)
    # print(img_hoechst.shape)
    img_green = img_crop(skio.imread("%s%s/%s/%s_CH2.tif" % (data_dir, folder, sample, file_name), plugin="tifffile")[:, :, 1], i)
    img_red = img_crop(skio.imread("%s%s/%s/%s_CH3.tif" % (data_dir, folder, sample, file_name), plugin="tifffile")[:, :, 0], i)
    img_farred = img_crop(skio.imread("%s%s/%s/%s_CH4.tif" % (data_dir, folder, sample, file_name), plugin="tifffile")[:, :, 0], i)

    if manual == 'YES':
        if os.path.exists("%s/%s/manual_seg/%s/%s_m-seg.tif" % (output_dir, folder, sample, file_name)):
            img_seg = skio.imread("%s%s/manual_seg/%s/%s_m-seg.tif" % (output_dir, folder, sample, file_name),
                                  plugin="tifffile")
        else:
            img_seg = np.zeros_like(img_hoechst)

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 30000])
        viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, 30000])
        viewer.add_image(img_red, blending='additive', colormap='red', contrast_limits=[0, 30000])
        viewer.add_image(img_farred, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
        viewer.add_image(img_seg, blending='additive', contrast_limits=[0, 1])
        shapes_add = viewer.add_shapes(name='add', ndim=2)
        napari.run()

        img_seg = erosion(img_seg, disk(6))
        img_seg = ima.napari_add_or_remove_obj(shapes_add.data, 'add', img_seg)
        img_seg = dilation(img_seg, disk(6))
        # img_seg = obj.label_remove_small_large_resort(label(img_seg), 169, 169)
        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 30000])
        viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, 30000])
        viewer.add_image(img_red, blending='additive', colormap='red', contrast_limits=[0, 30000])
        viewer.add_image(img_farred, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
        viewer.add_image(img_seg, blending='additive', contrast_limits=[0, 1])
        napari.run()
        if not os.path.exists("%s%s/manual_seg/%s/" % (output_dir, folder, sample)):
            os.makedirs("%s%s/manual_seg/%s/" % (output_dir, folder, sample))
        tif.imwrite("%s/%s/manual_seg/%s/%s_m-seg.tif" % (output_dir, folder, sample, file_name), img_seg)
    else:
        img_seg = skio.imread("%s%s/manual_seg/%s/%s_m-seg.tif" % (output_dir, folder, sample, file_name), plugin="tifffile")
    label_img = img_seg
    blue_props = regionprops(label_img, img_hoechst)
    green_props = regionprops(label_img, img_green)
    red_props = regionprops(label_img, img_red)
    farred_props = regionprops(label_img, img_farred)
    fov_lst = fov_lst + [i] * len(blue_props)
    hoechst_lst = hoechst_lst + [blue_props[j].intensity_mean for j in range(len(blue_props))]
    h23_lst = h23_lst + [green_props[j].intensity_mean for j in range(len(green_props))]
    azaleaB5_lst = azaleaB5_lst + [red_props[j].intensity_mean for j in range(len(red_props))]
    emiRFP670_lst = emiRFP670_lst + [farred_props[j].intensity_mean for j in range(len(farred_props))]

data['fov'] = fov_lst
data['hoechst'] = hoechst_lst
data['H2-3'] = h23_lst
data['AzaleaB5'] = azaleaB5_lst
data['emiRFP670'] = emiRFP670_lst
data['log10_emiRFP670'] = np.log10(data['emiRFP670'])

hc = 11000
cutoff = 3.1
n_total = len(data)
n_filtered = len(data[data['hoechst'] > hc])
n_neg = len(data[(data['hoechst'] > hc) & (data['log10_emiRFP670'] < cutoff)])
n_pos = len(data[(data['hoechst'] > hc) & (data['log10_emiRFP670'] >= cutoff)])
per_neg = n_neg/(n_filtered+0.01)
per_pos = n_pos/(n_filtered+0.01)

print(n_filtered)
print(n_neg)
print(n_pos)
print(per_neg)
print(per_pos)

data.to_csv('%s/%s/txt/%s_manual.txt' % (output_dir, folder, sample_name), index=False, sep='\t')
print("DONE!")
