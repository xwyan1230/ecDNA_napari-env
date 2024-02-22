import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
from skimage.measure import regionprops, label
from skimage.morphology import binary_dilation, disk, binary_erosion, dilation, erosion
import numpy as np
import shared.segmentation as seg
import shared.objects as obj
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240207_analysis_IFandRNAFISH_test/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

# set parameters
otsu_factor = 1.1
circ_threshold = 0.5
min_size = 300
max_size = 2000

folder = 'RNAFISH'
sample = 'mix_384'
sample_label = 'XY01'
n_img = 9

data = pd.DataFrame()
fov_lst = []
hoechst_lst = []
GFP_lst = []
mCherry_lst = []
RNAFISH_lst = []


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


for i in range(n_img):
    print(i)
    file_name = 'Image_%s_0000%s' % (sample_label, i+1)
    img_hoechst = img_crop(skio.imread("%s%s/%s/%s_CH1.tif" % (data_dir, folder, sample, file_name), plugin="tifffile")[:, :, 2], i)
    img_green = img_crop(skio.imread("%s%s/%s/%s_CH2.tif" % (data_dir, folder, sample, file_name), plugin="tifffile")[:, :, 1], i)
    img_red = img_crop(skio.imread("%s%s/%s/%s_CH3.tif" % (data_dir, folder, sample, file_name), plugin="tifffile")[:, :, 0], i)
    img_farred = img_crop(skio.imread("%s%s/%s/%s_CH4.tif" % (data_dir, folder, sample, file_name), plugin="tifffile")[:, :, 0], i)
    img_nuclear_seg = seg.cell_seg_fluorescent(img_hoechst, otsu_factor=otsu_factor, max_size=max_size, maxima_threshold=3,
                                               min_size=min_size, circ_thresh=circ_threshold).astype(int)
    print(np.max(img_nuclear_seg))
    img_nuclear_seg = erosion(img_nuclear_seg)
    nuclear_props = regionprops(img_nuclear_seg)
    RNA_mask = dilation(img_nuclear_seg, disk(5))
    RNA_mask[img_nuclear_seg >= 1] = 0

    blue_props = regionprops(img_nuclear_seg, img_hoechst)
    green_props = regionprops(img_nuclear_seg, img_green)
    red_props = regionprops(img_nuclear_seg, img_red)
    farred_props = regionprops(RNA_mask, img_farred)
    fov_lst = fov_lst + [i] * len(blue_props)
    hoechst_lst = hoechst_lst + [blue_props[j].intensity_mean for j in range(len(blue_props))]
    GFP_lst = GFP_lst + [green_props[j].intensity_mean for j in range(len(green_props))]
    mCherry_lst = mCherry_lst + [red_props[j].intensity_mean for j in range(len(red_props))]
    RNAFISH_lst = RNAFISH_lst + [farred_props[j].intensity_mean for j in range(len(farred_props))]
    print([red_props[j].label for j in range(len(red_props))])
    print([farred_props[j].label for j in range(len(farred_props))])

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(img_red, blending='additive', colormap='red', contrast_limits=[0, 65535])
    viewer.add_image(img_farred, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
    if not os.path.exists("%s%s/color_img/%s/" % (output_dir, folder, sample)):
        os.makedirs("%s%s/color_img/%s/" % (output_dir, folder, sample))
    plt.imsave("%s%s/color_img/%s/img_%s.tiff" % (output_dir, folder, sample, i+1), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_seg, blending='additive', colormap='red', contrast_limits=[0, np.max(img_nuclear_seg)])
    viewer.add_image(RNA_mask, blending='additive', colormap='green', contrast_limits=[0, np.max(RNA_mask)])
    plt.imsave("%s%s/color_img/%s/seg_%s.tiff" % (output_dir, folder, sample, i + 1), dis.blending(viewer))
    viewer.close()

data['fov'] = fov_lst
data['hoechst'] = hoechst_lst
data['GFP'] = GFP_lst
data['mCherry'] = mCherry_lst
data['RNAFISH'] = RNAFISH_lst
if not os.path.exists("%s%s/txt/" % (output_dir, folder)):
    os.makedirs("%s%s/txt/" % (output_dir, folder))
data.to_csv('%s/%s/txt/%s.txt' % (output_dir, folder, sample), index=False, sep='\t')
print("DONE!")