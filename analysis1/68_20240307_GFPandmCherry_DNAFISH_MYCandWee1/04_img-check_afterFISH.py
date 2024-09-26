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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

# set parameters
otsu_factor = 1.1
circ_threshold = 0.5
min_size = 300
max_size = 2000
# threshold = 3000

name = pd.read_csv('%s/name.txt' % data_dir, na_values=['.'], sep='\t')

k = 0
samples_lst = name['sample'].tolist()
wells_lst = name['well'].tolist()
sample = samples_lst[k]
well = wells_lst[k]
print(sample)
print(well)


img_hoechst = skio.imread("%s/afterFISH/%s_hoechst.tif" % (data_dir, well), plugin="tifffile")[:, :, 2]

viewer = napari.Viewer()
viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
# if not os.path.exists("%s%s/seg_update1/%s/" % (output_dir, folder, sample)):
#     os.makedirs("%s%s/seg_update1/%s/" % (output_dir, folder, sample))
# plt.imsave("%s%s/seg_update1/%s/seg_%s.tiff" % (output_dir, folder, sample, i+1), dis.blending(viewer))
# viewer.close()
napari.run()

print("DONE!")