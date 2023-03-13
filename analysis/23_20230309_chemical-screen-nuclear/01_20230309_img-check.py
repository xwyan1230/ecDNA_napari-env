import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230309_analysis_chemical-screen-nuclear/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

plate = 'HSR_2hr'
rows = ['B', 'C', 'D', 'E', 'F', 'G']
columns = ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
num_total_samples = len(rows) * len(columns)
for i in range(num_total_samples):
    row = rows[int(i/len(columns))]
    print(row)
    column = columns[int(i - int(i/len(columns))*len(columns))]
    print(column)
    imgs = [x for x in os.listdir('%s%s/%s/' % (data_dir, plate, row))]
    if '.DS_Store' in imgs:
        imgs.remove('.DS_Store')
    imgs = [x for x in imgs if '%s_%s%s' % (row, row, column) in x]
    n_imgs = int(len(imgs)/2)
    for j in range(n_imgs):
        file_name = '%s_%s%s_%s_RAW' % (row, row, column, j+1)
        img_hoechst = skio.imread("%s%s/%s/%s_ch01.tif" % (data_dir, plate, row, file_name), plugin="tifffile")
        img_MYC = skio.imread("%s%s/%s/%s_ch00.tif" % (data_dir, plate, row, file_name), plugin="tifffile")
        img_hoechst = np.concatenate([np.zeros(shape=[2048, 5]), img_hoechst], axis=1)[:2048, :2048]

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
        napari.run()
