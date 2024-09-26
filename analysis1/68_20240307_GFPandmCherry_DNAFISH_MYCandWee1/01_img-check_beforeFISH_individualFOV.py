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
n_img = 88

name = pd.read_csv('%s/name.txt' % data_dir, na_values=['.'], sep='\t')

samples_lst = name['sample'].tolist()
wells_lst = name['well'].tolist()
beforeFISH_lst = name['beforeFISH'].tolist()
afterFISH_lst = name['afterFISH'].tolist()

for k in range(len(samples_lst)):
    sample = samples_lst[k]
    well = wells_lst[k]
    beforeFISH = beforeFISH_lst[k]
    print(sample)
    print(well)

    for i in range(n_img):
        print(i)
        if i < 9:
            file_name = 'Image_%s_0000%s' % (beforeFISH, i+1)
        else:
            file_name = 'Image_%s_000%s' % (beforeFISH, i+1)
        img_green = skio.imread("%s/beforeFISH/%s/%s_CH2.tif" % (data_dir, beforeFISH, file_name), plugin="tifffile")[:, :, 1]
        img_red = skio.imread("%s/beforeFISH/%s/%s_CH3.tif" % (data_dir, beforeFISH, file_name), plugin="tifffile")[:, :, 0]

        viewer = napari.Viewer()
        viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, 65535])
        viewer.add_image(img_red, blending='additive', colormap='red', contrast_limits=[0, 65535])
        # if not os.path.exists("%s%s/seg_update1/%s/" % (output_dir, folder, sample)):
        #     os.makedirs("%s%s/seg_update1/%s/" % (output_dir, folder, sample))
        # plt.imsave("%s%s/seg_update1/%s/seg_%s.tiff" % (output_dir, folder, sample, i+1), dis.blending(viewer))
        # viewer.close()
        napari.run()

#     if not os.path.exists("%s%s/txt/" % (output_dir, folder)):
#         os.makedirs("%s%s/txt/" % (output_dir, folder))
#     data.to_csv('%s/%s/txt/%s_update1.txt' % (output_dir, folder, sample), index=False, sep='\t')
print("DONE!")