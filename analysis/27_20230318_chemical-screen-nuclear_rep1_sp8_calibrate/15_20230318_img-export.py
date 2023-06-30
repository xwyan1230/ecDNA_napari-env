import skimage.io as skio
import napari
from skimage.morphology import erosion, disk
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230318_analysis_chemical-screen-nuclear_rep1_sp8_calibrate/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

plate = 'HSR_6hr'
seq = pd.read_csv('%s%s_sequence.txt' % (data_dir, plate), na_values=['.'], sep='\t')

samples = seq['sample'].tolist()
start_sample = 54

n_nuclear_convex_dilation = -3

for s in range(len(samples)):
    print(s + start_sample)
    sample = samples[s + start_sample]
    print(sample)

    file_name = '20230318_HSR_6hr_rep1_calibration_%s_RAW' % (sample)
    img_hoechst_stack = skio.imread("%s%s/%s_ch01.tif" % (data_dir, plate, file_name), plugin="tifffile")
    img_MYC_stack = skio.imread("%s%s/%s_ch00.tif" % (data_dir, plate, file_name), plugin="tifffile")

    for fov in range(4):
        img_hoechst = img_hoechst_stack[:, :, fov]
        img_MYC = img_MYC_stack[:, :, fov]
        img_nuclear_seg = skio.imread("%s%s/seg_tif/%s_seg_%s.tif" % (data_dir1, plate, file_name, fov),
                                      plugin="tifffile")
        if n_nuclear_convex_dilation < 0:
            img_nuclear_seg = erosion(img_nuclear_seg, disk(abs(n_nuclear_convex_dilation)))
        # img_MYC = np.concatenate([np.zeros(shape=[5, 2048]), img_MYC], axis=0)[:2048, :2048]

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
        plt.imsave("%s%s/color_img/%s_img_%s.tiff" % (output_dir, plate, file_name, fov), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
        viewer.add_image(img_nuclear_seg, blending='additive')
        plt.imsave("%s%s/color_img/%s_seg_%s.tiff" % (output_dir, plate, file_name, fov), dis.blending(viewer))
        viewer.close()


