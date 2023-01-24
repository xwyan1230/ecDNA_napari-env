import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import random
import shared.image as ima
from skimage.morphology import dilation, disk
import napari
import tifffile as tif
import shared.display as dis
import os


def return_xy(r: int):
    out = []
    for i in range(r):
        for j in range(r):
            out.append((i, j))
    return out


def return_neighbour(dot, r):
    out = []
    neighbour = disk(r)
    for i in range(2*r+1):
        for j in range(2*r+1):
            if neighbour[i][j] == 1:
                if (i != (r+1)) | (j != (r+1)):
                    out.append((dot[0]+i-(r+1), dot[1]+j-(r+1)))
    return out


# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230119_model_clustering/"
data_dir = "%s" % master_folder
output_dir = "%s" % master_folder

r = 75
local_size = 150
coefficient = 100
krange = 5

im_test = np.zeros((2*local_size, 2*local_size))
im_seg = ima.logical_ellipse(im_test, local_size, local_size, r, r)
p_list = im_seg.flatten()
p_coordinate = return_xy(2*local_size)
p_dict = dict(zip(p_coordinate, p_list))

copy_num = random.choices(np.arange(20, 200, 1), k=1)
for k in range(len(copy_num)):
    print('coefficient: %s, cell: %s' % (coefficient, k))
    random_ones = []
    for i in range(copy_num[k]):
        weights = list(p_dict.values())
        dot = random.choices(p_coordinate, weights=weights, k=1)[0]
        random_ones.append(dot)
        p_dict[dot] = 0
        for j in return_neighbour(dot, krange):
            temp = p_dict[j]
            if temp != 0:
                p_dict[j] = temp + coefficient

    im_FISH = np.zeros((2 * local_size, 2 * local_size))
    for i in random_ones:
        im_FISH[i[0]][i[1]] = 1

    im_FISH = dilation(im_FISH, disk(2))
    im_FISH[im_seg == 0] = 0

    if not os.path.exists("%s%s_%s/seg_tif/" % (output_dir, coefficient, krange)):
        os.makedirs("%s%s_%s/seg_tif/" % (output_dir, coefficient, krange))
    tif.imwrite("%s%s_%s/seg_tif/k%s_r%s_%s_cp%s.tif" % (output_dir, coefficient, krange, coefficient, krange, k, copy_num[k]), im_FISH)

    viewer = napari.Viewer()
    viewer.add_image(im_seg, colormap='blue', blending='additive')
    viewer.add_image(im_FISH, colormap='green', blending='additive', contrast_limits=[0, 1])
    if not os.path.exists("%s%s_%s/color_tif/" % (output_dir, coefficient, krange)):
        os.makedirs("%s%s_%s/color_tif/" % (output_dir, coefficient, krange))
    plt.imsave('%s%s_%s/color_tif/k%s_r%s_%s_cp%s.tiff' % (output_dir, coefficient, krange, coefficient, krange, k, copy_num[k]), dis.blending(viewer))
    viewer.close()

print("DONE!")