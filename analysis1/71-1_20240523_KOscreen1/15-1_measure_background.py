import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import imutils
import pandas as pd
import shared.image as ima
import tifffile as tif
import math
import nd2
import shared.segmentation as seg
import shared.objects as obj
import os
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'G2'
total_fov = 16
dshape_factor = 0.145
pixel_size = 300/2720  # uM
img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, sample))

data = pd.DataFrame(columns=['sample', 'fov', 'bg_merge', 'bg_z'])
pd_seg_z = pd.read_csv('%s/%s/10_%s_seg_z.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')

for fov in range(total_fov):
    print("%s/%s" % (fov+1, total_fov))
    img_hoechst = img_stack[fov, :, 0, :, :]
    img_hoechst_merge = img_hoechst.max(axis=0)
    img_DNAFISH = img_stack[fov, :, 1, :, :]
    img_DNAFISH_sum = img_DNAFISH.astype(int).sum(axis=0)
    seg_z = pd_seg_z['seg_z'][fov]
    img_DNAFISH_seg_z = img_DNAFISH[seg_z]

    img_seg_bg = np.zeros_like(img_hoechst_merge)

    viewer = napari.Viewer()
    # viewer.add_image(img_hoechst_merge, blending='additive', colormap='blue', contrast_limits=[0, 12000])
    viewer.add_image(img_DNAFISH_sum, blending='additive', colormap='green', contrast_limits=[0, 65535])
    shapes = viewer.add_shapes(name='Shapes', ndim=2)
    napari.run()

    img_seg_bg = ima.napari_add_or_remove_obj(shapes.data, 'add', img_seg_bg)

    if not os.path.exists("%s/%s/15_seg_out_bg_tif/" % (output_dir, sample)):
        os.makedirs("%s/%s/15_seg_out_bg_tif/" % (output_dir, sample))
    tif.imwrite("%s/%s/15_seg_out_bg_tif/%s_%s_seg_out_bg.tif" % (output_dir, sample, sample, fov), img_seg_bg)

    total_props = regionprops(label(img_seg_bg), img_DNAFISH_sum)
    z_props = regionprops(label(img_seg_bg), img_DNAFISH_seg_z)

    data_temp = pd.DataFrame()
    data_temp['sample'] = [sample] * len(total_props)
    data_temp['fov'] = [fov] * len(total_props)
    data_temp['bg_merge'] = [total_props[i].intensity_mean for i in range(len(total_props))]
    data_temp['bg_z'] = [z_props[i].intensity_mean for i in range(len(z_props))]
    data = pd.concat([data, data_temp], axis=0)

data.to_csv('%s/%s/15_%s_bg.txt' % (output_dir, sample, sample), index=False, sep='\t')

if os.path.exists("%s/bg.txt" % output_dir):
    df_bg = pd.read_csv('%s/bg.txt' % output_dir, na_values=['.'], sep='\t')
else:
    df_bg = pd.DataFrame(columns=['sample', 'bg_merge', 'bg_z', 'n_z'])
if len(df_bg[df_bg['sample'] == sample]) == 0:
    df_bg.loc[len(df_bg.index)] = [sample, np.mean(data['bg_merge']), np.mean(data['bg_z']), img_stack.shape[1]]
else:
    location = df_bg[df_bg['sample'] == sample].index
    df_bg.loc[location] = [sample, np.mean(data['bg_merge']), np.mean(data['bg_z']), img_stack.shape[1]]
df_bg.to_csv('%s/bg.txt' % output_dir, index=False, sep='\t')

print("DONE!")
