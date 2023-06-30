import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import shared.objects as obj
import tifffile as tif
from skimage.morphology import erosion, disk
from skimage.measure import label, regionprops
import math
import shared.segmentation as seg
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230313_analysis_chemical-screen-nuclear_rep1_EVOS/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

# set parameters
otsu_factor = 1.2
circ_threshold = 0.8
min_size = 400
max_size = 2000

plate = 'DM_6hr'
rows = ['G']
columns = ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
num_total_samples = len(rows) * len(columns)
for i in range(num_total_samples):
    row = rows[int(i/len(columns))]
    column = columns[int(i - int(i/len(columns))*len(columns))]
    print('%s%s' % (row, column))
    imgs = [x for x in os.listdir('%s%s/' % (data_dir, plate))]
    if '.DS_Store' in imgs:
        imgs.remove('.DS_Store')
    imgs = [x for x in imgs if '%s%s' % (row, column) in x]
    n_imgs = int(len(imgs)/3)
    for j in range(n_imgs):
        file_name = '%s%s _000%s' % (row, column, j+1)
        print(file_name)
        img_hoechst = skio.imread("%s%s/%s_tgBFP.tif" % (data_dir, plate, file_name), plugin="tifffile")  # 1536, 2048
        img_MYC = skio.imread("%s%s/%s_CY5.tif" % (data_dir, plate, file_name), plugin="tifffile")

        # nuclear seg
        img_nuclear_seg = seg.cell_seg_fluorescent(img_hoechst, otsu_factor=otsu_factor, max_size=max_size, min_size=min_size, circ_thresh=circ_threshold).astype(int)

        if not os.path.exists("%s%s/seg_tif/" % (output_dir, plate)):
            os.makedirs("%s%s/seg_tif/" % (output_dir, plate))
        tif.imwrite("%s%s/seg_tif/%s_seg.tif" % (output_dir, plate, file_name), img_nuclear_seg)

        if not os.path.exists("%s%s/color_img/" % (output_dir, plate)):
            os.makedirs("%s%s/color_img/" % (output_dir, plate))

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 255])
        viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 255])
        viewer.add_image(img_nuclear_seg, blending='additive')
        plt.imsave("%s%s/color_img/%s_seg.tiff" % (output_dir, plate, file_name), dis.blending(viewer))
        viewer.close()
        # napari.run()
