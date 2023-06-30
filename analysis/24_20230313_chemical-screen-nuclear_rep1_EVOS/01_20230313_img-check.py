import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230313_analysis_chemical-screen-nuclear_rep1_EVOS/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

plate = 'HSR_6hr'
rows = ['B', 'C', 'D', 'E', 'F', 'G']
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
        img_hoechst = skio.imread("%s%s/%s_tgBFP.tif" % (data_dir, plate, file_name), plugin="tifffile")  # 1536, 2048
        img_MYC = skio.imread("%s%s/%s_CY5.tif" % (data_dir, plate, file_name), plugin="tifffile")
        # img_hoechst = np.concatenate([np.zeros(shape=[2048, 5]), img_hoechst], axis=1)[:2048, :2048]

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 255])
        viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 255])
        napari.run()
