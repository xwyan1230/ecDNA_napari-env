import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230320_analysis_chemical-screen-nuclear_DNAFISH/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

plate = 'DM_24hr'
samples = ['D9', 'C3']
"""seq = pd.read_csv('%s%s_sequence.txt' % (data_dir, plate), na_values=['.'], sep='\t')

samples = seq['sample'].tolist()
start_sample = 27"""

"""imgs = [x for x in os.listdir('%s%s/' % (data_dir, plate))]
if '.DS_Store' in imgs:
    imgs.remove('.DS_Store')
imgs1 = [imgs[i].split('_RAW')[0].split('calibrate_')[1] for i in range(len(imgs))]
for i in imgs1:
    if i not in samples:
        print(i)"""

for s in range(len(samples)):
    sample = samples[s]

    # file_name = '20230318_HSR_6hr_rep1_calibration_%s_RAW' % (sample)
    file_name = '20230317_DM_24hr_DNAFISH_%s_%s_RAW' % (sample[0], sample)
    img_hoechst_stack = skio.imread("%s%s/%s_ch01.tif" % (data_dir, plate, file_name), plugin="tifffile")
    img_DNAFISH_stack = skio.imread("%s%s/%s_ch00.tif" % (data_dir, plate, file_name), plugin="tifffile")

    for fov in range(9):
        img_hoechst = img_hoechst_stack[fov, :, :]
        img_DNAFISH = img_DNAFISH_stack[fov, :, :]
        # img_MYC = np.concatenate([np.zeros(shape=[5, 2048]), img_MYC], axis=0)[:2048, :2048]

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 45535])
        viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
        plt.imsave("%s%s/%s_img_%s.tiff" % (output_dir, plate, file_name, fov), dis.blending(viewer))
        napari.run()


