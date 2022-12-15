import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221214_analysis_Natasha_GBMec_dual-funnel/221017_GBMEC_3hrtreatment_EGFR_interphaseFISH/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'slide1_DMSO'

start_fov = 1
total_fov = 50

for f in range(total_fov):
    fov = f + start_fov
    file_name = 'TileScan 1_Position %s_RAW' % fov
    img_hoechst = skio.imread("%s221017_GBMEC_%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_DNAFISH = skio.imread("%s221017_GBMEC_%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    napari.run()
