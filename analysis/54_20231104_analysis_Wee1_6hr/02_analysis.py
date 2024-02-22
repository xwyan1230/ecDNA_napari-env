import skimage.io as skio
import napari
import os
import tifffile as tif
import shared.image as ima
from skimage.measure import label, regionprops
from skimage.morphology import binary_dilation, disk
import numpy as np
import pandas as pd

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231104_analysis_ColoDMandHSR_Wee1_6hr/"

sample = 'ColoHSR_Wee1_degrader_gH2AX_pH3'
conc = 'XY42'

# img
img_CH1 = skio.imread("%s/data/%s/%s/Image_%s_CH1.tif" % (master_folder, sample, conc, conc),
                      plugin="tifffile")[:, :, 2]
img_CH2 = skio.imread("%s/data/%s/%s/Image_%s_CH2.tif" % (master_folder, sample, conc, conc),
                      plugin="tifffile")[:, :, 1]
img_CH3 = skio.imread("%s/data/%s/%s/Image_%s_CH3.tif" % (master_folder, sample, conc, conc),
                      plugin="tifffile")[:, :, 0]
img_CH4 = skio.imread("%s/data/%s/%s/Image_%s_CH4.tif" % (master_folder, sample, conc, conc),
                      plugin="tifffile")[:, :, 0]
img_seg = skio.imread("%s/data/%s/%s/Image_%s_seg.tif" % (master_folder, sample, conc, conc),
                      plugin="tifffile")

viewer = napari.Viewer()
viewer.add_image(img_CH1, blending='additive', colormap='blue', contrast_limits=[0, 30000])
viewer.add_image(img_CH2, blending='additive', colormap='green', contrast_limits=[0, 30000])
viewer.add_image(img_CH3, blending='additive', colormap='red', contrast_limits=[0, 30000])
viewer.add_image(img_CH4, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
viewer.add_image(img_seg, blending='additive', contrast_limits=[0, 1])
napari.run()

gH2AX = regionprops(label(img_seg), img_CH1)
green = regionprops(label(img_seg), img_CH2)
red = regionprops(label(img_seg), img_CH3)
pH3 = regionprops(label(img_seg), img_CH4)

data = pd.DataFrame()
data['gH2AX'] = [gH2AX[i].intensity_mean for i in range(len(gH2AX))]
data['green'] = [green[i].intensity_mean for i in range(len(green))]
data['red'] = [red[i].intensity_mean for i in range(len(red))]
data['pH3'] = [pH3[i].intensity_mean for i in range(len(pH3))]

data.to_csv('%s/%s_%s.txt' % (master_folder, sample, conc), index=False, sep='\t')