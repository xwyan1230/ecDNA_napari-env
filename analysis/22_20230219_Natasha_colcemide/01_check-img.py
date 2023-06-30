import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import tifffile as tif
import shared.image as ima

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230209_analysis_Natasha_DMcolcemid_mchctrl/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

start_fov = 1
end_fov = 10
seq = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10']

sample = 'DMcolcemid_mchctrl'
total_fov = end_fov - start_fov + 1

for f in range(total_fov):
    fov = start_fov + f
    print(fov)
    file_name = '230209_%s_%s-OME TIFF-Export-%s' % (sample, fov, seq[f])
    img_hoechst = skio.imread("%s%s/%s.ome.tiff" % (data_dir, sample, file_name), plugin="tifffile")[:, :, 0]
    img_DNAFISH = skio.imread("%s%s/%s.ome.tiff" % (data_dir, sample, file_name), plugin="tifffile")[:, :, 1]
    img_mCherry = skio.imread("%s%s/%s.ome.tiff" % (data_dir, sample, file_name), plugin="tifffile")[:, :, 2]
    print(img_hoechst.shape)

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    viewer.add_image(img_mCherry, blending='additive', colormap='red', contrast_limits=[0, img_mCherry.max()])
    napari.run()






