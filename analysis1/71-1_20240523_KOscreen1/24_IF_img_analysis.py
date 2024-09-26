import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
from skimage.measure import regionprops, label
from skimage.morphology import binary_dilation, disk, binary_erosion, dilation, erosion
import numpy as np
import vispy.color
import shared.segmentation as seg
import shared.objects as obj
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

# set parameters
otsu_factor = 1.1
circ_threshold = 0.5
min_size = 300
max_size = 2000
threshold = 3000

folder = 'IF'
sample = 'G2'
rep = 2
df = pd.read_csv('%s/koscreen1.txt' % master_folder, na_values=['.'], sep='\t')
ctrls = ['J3', 'J4', 'J5', 'J6', 'J7', 'J8', 'J9', 'J10', 'J11', 'J12']
if sample not in ctrls:
    sample_label = df[df['96well'] == sample]['384well_rep%s_%s' % (rep, folder)].tolist()[0]
else:
    sample_label = df[df['384well_rep%s_well' % rep] == sample]['384well_rep%s_%s' % (rep, folder)].tolist()[0]
print(sample_label)
n_img = 9

data = pd.DataFrame()
fov_lst = []
hoechst_lst = []
GFP_lst = []
mCherry_lst = []
IF_lst = []


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
    img_hoechst = img_crop(skio.imread("%s%s/rep%s/%s/%s_CH1.tif" % (data_dir, folder, rep, sample_label, file_name), plugin="tifffile")[:, :, 2], i)
    img_green = img_crop(skio.imread("%s%s/rep%s/%s/%s_CH2.tif" % (data_dir, folder, rep, sample_label, file_name), plugin="tifffile")[:, :, 1], i)
    img_red = img_crop(skio.imread("%s%s/rep%s/%s/%s_CH3.tif" % (data_dir, folder, rep, sample_label, file_name), plugin="tifffile")[:, :, 0], i)
    img_farred = img_crop(skio.imread("%s%s/rep%s/%s/%s_CH4.tif" % (data_dir, folder, rep, sample_label, file_name), plugin="tifffile")[:, :, 0], i)

    img_nuclear_seg = seg.cell_seg_fluorescent(img_hoechst, otsu_factor=otsu_factor, max_size=max_size,
                                               maxima_threshold=1.0001,
                                               min_size=min_size, circ_thresh=circ_threshold,
                                               threshold=threshold).astype(int)
    print(np.max(img_nuclear_seg))
    img_nuclear_seg = erosion(img_nuclear_seg, disk(4))
    nuclear_props = regionprops(img_nuclear_seg)
    centroids = [nuclear_props[j].centroid for j in range(len(nuclear_props))]
    img_seg = np.zeros_like(img_nuclear_seg)
    for j in range(len(centroids)):
        img_seg[int(centroids[j][0])][int(centroids[j][1])] = 1
    img_seg = binary_dilation(img_seg, disk(7))
    img_seg = obj.label_remove_small_large_resort(label(img_seg), 149, 149)
    label_img = label(img_seg)
    # RNA_mask = dilation(img_nuclear_seg, disk(5))
    # RNA_mask[img_nuclear_seg >= 1] = 0

    blue_props = regionprops(label_img, img_hoechst)
    green_props = regionprops(label_img, img_green)
    red_props = regionprops(label_img, img_red)
    farred_props = regionprops(label_img, img_farred)
    fov_lst = fov_lst + [i] * len(blue_props)
    hoechst_lst = hoechst_lst + [blue_props[j].intensity_mean for j in range(len(blue_props))]
    GFP_lst = GFP_lst + [green_props[j].intensity_mean for j in range(len(green_props))]
    mCherry_lst = mCherry_lst + [red_props[j].intensity_mean for j in range(len(red_props))]
    IF_lst = IF_lst + [farred_props[j].intensity_mean for j in range(len(farred_props))]

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(img_red, blending='additive', colormap='red', contrast_limits=[0, 65535])
    viewer.add_image(img_farred, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
    if not os.path.exists("%s%s/24_IF_rep%s_color_img/" % (output_dir, sample, rep)):
        os.makedirs("%s%s/24_IF_rep%s_color_img/" % (output_dir, sample, rep))
    plt.imsave("%s%s/24_IF_rep%s_color_img/%s_IF_img_%s.tiff" % (output_dir, sample, rep, sample, i+1), dis.blending(viewer))
    viewer.close()

    n_nuclear = len(nuclear_props)
    cmap = plt.cm.get_cmap('viridis')
    colors = []
    for k in range(n_nuclear):
        colors.append(cmap((k+1)/n_nuclear)[:3])
    np.random.shuffle(colors)
    random_colors = [[0, 0, 0]] + colors

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(label_img, blending='additive', colormap=vispy.color.Colormap(random_colors),
                     contrast_limits=[0, np.max(label_img)])
    # viewer.add_image(img_nuclear_seg, blending='additive', colormap=vispy.color.Colormap(random_colors), contrast_limits=[0, np.max(img_nuclear_seg)])
    plt.imsave("%s%s/24_IF_rep%s_color_img/%s_seg_%s.tiff" % (output_dir, sample, rep, sample, i + 1), dis.blending(viewer))
    viewer.close()

data['fov'] = fov_lst
data['hoechst'] = hoechst_lst
data['GFP'] = GFP_lst
data['mCherry'] = mCherry_lst
data['IF'] = IF_lst
data.to_csv('%s/%s/24_%s_IF_rep%s.txt' % (output_dir, sample, sample, rep), index=False, sep='\t')
print("DONE!")