import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
from skimage.measure import regionprops, label
from skimage.morphology import binary_dilation, disk
import numpy as np
import shared.segmentation as seg
import shared.objects as obj
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240709_analysis_CRISPR-C/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

# set parameters
otsu_factor = 1.05
circ_threshold = 0.5
min_size = 300
max_size = 3000
threshold = 0


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


def img_factor(i):
    if i == 0:
        out = 1
    elif (i == 1) | (i == 2) | (i == 3) | (i == 6):
        out = 0.7
    elif (i == 4) | (i == 5) | (i == 7) | (i == 8):
        out = 0.49
    return out


sample = 'XY08'

n_img = 9

data = pd.DataFrame()
fov_lst = []
hoechst_lst = []
GFP_lst = []
fov_hoechst_lst = []
fov_GFP_lst = []

for i in range(n_img):
    print(i)
    file_name = 'Image_%s_0000%s' % (sample, i+1)
    img_hoechst = img_crop(skio.imread("%s%s/%s_CH1.tif" % (data_dir, sample, file_name), plugin="tifffile")[:, :, 2], i)
    img_green = img_crop(skio.imread("%s%s/%s_CH2.tif" % (data_dir, sample, file_name), plugin="tifffile")[:, :, 1], i)
    img_nuclear_seg = seg.cell_seg_fluorescent(img_hoechst, otsu_factor=otsu_factor, max_size=max_size,
                                               maxima_threshold=1.1,
                                               min_size=min_size, circ_thresh=circ_threshold,
                                               threshold=threshold).astype(int)
    fov_hoechst = np.sum(img_hoechst) / img_factor(i)
    fov_GFP = np.sum(img_green)/ img_factor(i)
    nuclear_props = regionprops(img_nuclear_seg)
    centroids = [nuclear_props[j].centroid for j in range(len(nuclear_props))]
    img_seg = np.zeros_like(img_nuclear_seg)
    for j in range(len(centroids)):
        img_seg[int(centroids[j][0])][int(centroids[j][1])] = 1
    img_seg = binary_dilation(img_seg, disk(7))
    img_seg = obj.label_remove_small_large_resort(label(img_seg), 149, 149)
    label_img = label(img_seg)
    blue_props = regionprops(label_img, img_hoechst)
    green_props = regionprops(label_img, img_green)

    fov_lst = fov_lst + [i] * len(blue_props)
    hoechst_lst = hoechst_lst + [blue_props[j].intensity_mean for j in range(len(blue_props))]
    GFP_lst = GFP_lst + [green_props[j].intensity_mean for j in range(len(green_props))]
    fov_hoechst_lst = fov_hoechst_lst + [fov_hoechst]* len(blue_props)
    fov_GFP_lst = fov_GFP_lst + [fov_GFP] * len(blue_props)

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(img_seg, blending='additive', contrast_limits=[0, 1])
    if not os.path.exists("%s/seg/%s/" % (output_dir, sample)):
        os.makedirs("%s/seg/%s/" % (output_dir, sample))
    plt.imsave("%s/seg/%s/seg_%s.tiff" % (output_dir, sample, i + 1), dis.blending(viewer))
    viewer.close()
    # napari.run()
data['sample'] = [sample] * len(fov_lst)
data['fov'] = fov_lst
data['hoechst'] = hoechst_lst
data['GFP'] = GFP_lst
data['fov_hoechst'] = fov_hoechst_lst
data['fov_GFP'] = fov_GFP_lst
if not os.path.exists("%s/txt/" % (output_dir)):
    os.makedirs("%s/txt/" % (output_dir))
data.to_csv('%s/txt/%s.txt' % (output_dir, sample), index=False, sep='\t')
print("DONE!")