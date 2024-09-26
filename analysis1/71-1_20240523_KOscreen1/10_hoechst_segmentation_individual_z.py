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

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'G2'
total_fov = 16
convex_conversion_threshold = 0.8
circ_threshold = 0.8
local_size = 200
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = 3000
max_size_nuclear = 9000
img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, sample))

data = pd.DataFrame(columns=['sample', 'fov', 'seg_z'])

for fov in range(total_fov):
    print("%s/%s" % (fov+1, total_fov))
    img_hoechst = img_stack[fov, :, 0, :, :]
    total_hoechst_lst = []
    for z in range(img_hoechst.shape[0]):
        total_hoechst = np.sum(img_hoechst[z])
        total_hoechst_lst.append(total_hoechst)
    seg_z = total_hoechst_lst.index(max(total_hoechst_lst))
    print(seg_z)
    data.loc[len(data.index)] = [sample, fov, seg_z]

    img_hoechst_seg = img_hoechst[seg_z]
    img_nuclear_seg = seg.nuclear_seg_nikon(img_hoechst_seg, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                            max_size=max_size_nuclear)
    img_nuclear_seg_convex = obj.label_remove_low_circ(seg.obj_to_convex_filter(img_nuclear_seg, threshold=convex_conversion_threshold), thresh=circ_threshold)

    if not os.path.exists("%s/%s/10_seg_z_tif/" % (output_dir, sample)):
        os.makedirs("%s/%s/10_seg_z_tif/" % (output_dir, sample))
    if not os.path.exists("%s/%s/10_seg_z_color/" % (output_dir, sample)):
        os.makedirs("%s/%s/10_seg_z_color/" % (output_dir, sample))
    tif.imwrite("%s/%s/10_seg_z_tif/%s_%s_seg_z.tif" % (output_dir, sample, sample, fov), img_nuclear_seg_convex)

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst_seg, blending='additive', colormap='blue', contrast_limits=[0, 12000])
    # viewer.add_image(img_nuclear_seg, blending='additive', colormap='green', contrast_limits=[0, 1])
    viewer.add_image(img_nuclear_seg_convex, blending='additive', contrast_limits=[0, 1])
    plt.imsave("%s%s/10_seg_z_color/%s_%s_seg_z_color.tiff" % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

data.to_csv('%s/%s/10_%s_seg_z.txt' % (output_dir, sample, sample), index=False, sep='\t')

print("DONE!")