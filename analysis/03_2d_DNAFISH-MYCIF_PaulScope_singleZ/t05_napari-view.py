import napari
import skimage.io as skio
import matplotlib.pyplot as plt
import os
import numpy as np

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220526_flowFISH_topHits_screen/"
sample = 'E9'
raw_folder = '01_raw'
save_folder = 'v5_raw-image'
save_path = '%s%s/' % (master_folder, save_folder)
if not os.path.exists(save_path):
    os.makedirs(save_path)

fov = 1
z = 17
blue_thresh = 65535
green_thresh = 6000
far_thresh = 10000


def blending(view):
    blended = np.zeros(view.layers[0].data.shape + (4,))
    for layer in view.layers:
        # normalize data by clims
        normalized_data = (layer.data - layer.contrast_limits[0]) / (layer.contrast_limits[1] - layer.contrast_limits[0])
        colormapped_data = layer.colormap.map(normalized_data.flatten())
        colormapped_data = colormapped_data.reshape(normalized_data.shape + (4,))

        blended = blended + colormapped_data

    blended[..., 3] = 1  # set alpha channel to 1
    blended[blended > 1] = 1
    blended[blended < 0] = 0
    return blended


im_z_stack_nuclear = skio.imread("%s%s/%s/%s/R%s_RAW_ch00.tif" %
                                 (master_folder, sample[0], sample[1:], raw_folder, fov), plugin="tifffile")
im_z_stack_DNAFISH = skio.imread("%s%s/%s/%s/R%s_RAW_ch01.tif" %
                                 (master_folder, sample[0], sample[1:], raw_folder, fov), plugin="tifffile")
im_z_stack_IF = skio.imread("%s%s/%s/%s/R%s_RAW_ch02.tif" %
                            (master_folder, sample[0], sample[1:], raw_folder, fov), plugin="tifffile")

"""viewer = napari.view_image(im_z_stack_nuclear, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(im_z_stack_DNAFISH, blending='additive', colormap='green', gamma=1.64, contrast_limits=[500, 8000])
viewer.add_image(im_z_stack_IF, blending='additive', colormap='magenta', gamma=1.68, contrast_limits=[0, 20000])
napari.run()"""

viewer3 = napari.view_image(im_z_stack_nuclear[z], blending='additive', colormap='blue', contrast_limits=[0, blue_thresh])
viewer3.add_image(im_z_stack_DNAFISH[z], blending='additive', colormap='green', contrast_limits=[0, green_thresh])
viewer3.add_image(im_z_stack_IF[z], blending='additive', colormap='magenta', contrast_limits=[0, far_thresh])
plt.imsave('%s%s/%s_merge_fov%s_z%s_b%s_g%s_f%s.tiff' % (master_folder, save_folder, sample, fov, z, blue_thresh,
                                                         green_thresh, far_thresh), blending(viewer3))

viewer = napari.view_image(im_z_stack_nuclear[z], blending='additive', colormap='blue', contrast_limits=[0, blue_thresh])
plt.imsave('%s%s/%s_nuclear_fov%s_z%s_b%s.tiff' % (master_folder, save_folder, sample, fov, z, blue_thresh), blending(viewer))
viewer1 = napari.view_image(im_z_stack_DNAFISH[z], blending='additive', colormap='green', contrast_limits=[0, green_thresh])
plt.imsave('%s%s/%s_DNAFISH_fov%s_z%s_g%s.tiff' % (master_folder, save_folder, sample, fov, z, green_thresh), blending(viewer1))
viewer2 = napari.view_image(im_z_stack_IF[z], blending='additive', colormap='magenta', contrast_limits=[0, far_thresh])
plt.imsave('%s%s/%s_IF_fov%s_z%s_f%s.tiff' % (master_folder, save_folder, sample, fov, z, far_thresh), blending(viewer2))


print("DONE!")