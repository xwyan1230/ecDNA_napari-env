import skimage.io as skio
import napari
import imutils
import shared.image as ima
import matplotlib.pyplot as plt
import numpy as np
import shared.display as dis
# from napari_animation import Animation
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230711_analysis_H3K4me3_NPC/"
data_dir = "%sdata/raw/" % master_folder
output_dir = "%sdata/" % master_folder

# sample = 'Colo320DM'
# sample = 'NPC'
sample = 'TPR'
# prefix = '20230616_H3K4me3_DM_H3K4me3'
# prefix = '20230608_DM_IF-NPC_acidFISH-MYC_ColoDM_IF-NPC_acidFISH-MYC_100x-2048'
prefix = '20230616_DM_TPR-Nup98_DM_TPR'
total_fov = 2
start_fov = 1
fov = 1
center = [1543.9267693 ,  904.96859135]
local_range = 200
position = [int(center[0]-local_range), int(center[0]+local_range), int(center[1]-local_range), int(center[1]+local_range)]

# The microscope reports the following spacing (in µm)
original_spacing = np.array([0.3, 0.05679, 0.05679])
# We downsampled each slice 4x to make the data smaller
rescaled_spacing = original_spacing * [1, 4, 4]
# Normalize the spacing so that pixels are a distance of 1 apart
spacing = rescaled_spacing / rescaled_spacing[2]

lst = [x.split('%s_' % prefix)[1].split('_ch')[0] for x in os.listdir('%s%s' % (data_dir, sample)) if x[-4:] == '.tif']

lst_temp = [int(x.split('_z')[1]) for x in lst if x.split('_z')[0] == '%s' % fov]
total_z = max(lst_temp) + 1
img_3d_nuclear = skio.imread("%s%s/%s_%s_z00_ch00.tif" % (data_dir, sample, prefix, fov), plugin="tifffile")[np.newaxis, position[0]:position[1], position[2]:position[3]]
img_3d_DNAFISH = skio.imread("%s%s/%s_%s_z00_ch01.tif" % (data_dir, sample, prefix, fov), plugin="tifffile")[np.newaxis, position[0]:position[1], position[2]:position[3]]
img_3d_laminB = skio.imread("%s%s/%s_%s_z00_ch02.tif" % (data_dir, sample, prefix, fov), plugin="tifffile")[np.newaxis, position[0]:position[1], position[2]:position[3]]
for z in range(total_z-1):
    if z < 9:
        file_name = '%s_%s_z0%s' % (prefix, fov, z+1)
    else:
        file_name = '%s_%s_z%s' % (prefix, fov, z+1)
    img_nuclear = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_DNAFISH = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_laminB = skio.imread("%s%s/%s_ch02.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_3d_nuclear = np.vstack([img_3d_nuclear, img_nuclear[np.newaxis, position[0]:position[1], position[2]:position[3]]])
    img_3d_DNAFISH = np.vstack([img_3d_DNAFISH, img_DNAFISH[np.newaxis, position[0]:position[1], position[2]:position[3]]])
    img_3d_laminB = np.vstack([img_3d_laminB, img_laminB[np.newaxis, position[0]:position[1], position[2]:position[3]]])
    print(img_3d_nuclear.shape)
viewer = napari.Viewer()
viewer.add_image(img_3d_nuclear, blending='additive', colormap='blue', name='nuclei', scale=original_spacing, contrast_limits=[0, 65535])
viewer.add_image(img_3d_DNAFISH, blending='additive', colormap='green', name='DNAFISH', scale=original_spacing, contrast_limits=[0, 45000])
viewer.add_image(img_3d_laminB, blending='additive', colormap='red', name='laminB', scale=original_spacing, contrast_limits=[0, 65535])
viewer.scale_bar.visible = True
viewer.scale_bar.unit = "um"
# viewer.window.add_plugin_dock_widget(plugin_name='napari-animation')
napari.run()

if not os.path.exists('%s/color_img/%s/fov%s_%s_%s/' % (output_dir, sample, fov, center[0], center[1])):
    os.makedirs('%s/color_img/%s/fov%s_%s_%s/' % (output_dir, sample, fov, center[0], center[1]))

for z in range(total_z):
    viewer = napari.Viewer()
    viewer.add_image(img_3d_nuclear[z, :, :], blending='additive', colormap='blue', name='nuclei', contrast_limits=[0, 65535])
    viewer.add_image(img_3d_DNAFISH[z, :, :], blending='additive', colormap='green', name='DNAFISH', contrast_limits=[0, 45000])
    viewer.add_image(img_3d_laminB[z, :, :], blending='additive', colormap='red', name='laminB', contrast_limits=[0, 65535])
    plt.imsave('%s/color_img/%s/fov%s_%s_%s/z%s_overlay.tiff' % (output_dir, sample, fov, center[0], center[1], z), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_3d_nuclear[z, :, :], blending='additive', colormap='blue', name='nuclei',contrast_limits=[0, 65535])
    viewer.add_image(img_3d_DNAFISH[z, :, :], blending='additive', colormap='green', name='DNAFISH',contrast_limits=[0, 45000])
    plt.imsave('%s/color_img/%s/fov%s_%s_%s/z%s_woLamin.tiff' % (output_dir, sample, fov, center[0], center[1], z), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_3d_nuclear[z, :, :], blending='additive', colormap='blue', name='nuclei', contrast_limits=[0, 65535])
    plt.imsave('%s/color_img/%s/fov%s_%s_%s/z%s_nuclear.tiff' % (output_dir, sample, fov, center[0], center[1], z), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_3d_DNAFISH[z, :, :], blending='additive', colormap='green', name='DNAFISH', contrast_limits=[0, 45000])
    plt.imsave('%s/color_img/%s/fov%s_%s_%s/z%s_DNAFISH.tiff' % (output_dir, sample, fov, center[0], center[1], z), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_3d_laminB[z, :, :], blending='additive', colormap='red', name='laminB', contrast_limits=[0, 65535])
    plt.imsave('%s/color_img/%s/fov%s_%s_%s/z%s_laminB.tiff' % (output_dir, sample, fov, center[0], center[1], z), dis.blending(viewer))
    viewer.close()

print("DONE!")

