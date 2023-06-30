import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230216_analysis_BRD4-series/"
data_dir = "%sdata/IF/20230210_KO-series_MYC-IF/" % master_folder
output_dir = "%sfigures/" % master_folder
start_fov = 1
end_fov = 6

sample = 'DM-CTCF-KO_oldMYC'
total_fov = end_fov - start_fov + 1

for f in range(total_fov):
    fov = start_fov + f
    file_name = '%s_Image %s_RAW' % (sample, fov)
    img_hoechst = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_MYC = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_MYC = np.concatenate([np.zeros(shape=[10, 3144]), img_MYC], axis=0)

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
    napari.run()
