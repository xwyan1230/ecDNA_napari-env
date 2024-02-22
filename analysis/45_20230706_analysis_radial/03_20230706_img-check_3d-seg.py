import skimage.io as skio
import napari
import imutils
import shared.image as ima
import numpy as np
import tifffile as tif
# from napari_animation import Animation
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230706_analysis_radial/"
data_dir = "%sdata/raw/" % master_folder
data_dir1 = "%sdata/seg_tif/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'Colo320DM_acidFISH_lamin_3d'
prefix = '20230614_acidFISH_lamin_ColoDM_DM'
total_fov = 8
start_fov = 1

# The microscope reports the following spacing (in µm)
original_spacing = np.array([0.3, 0.05679, 0.05679])
# We downsampled each slice 4x to make the data smaller
rescaled_spacing = original_spacing * [1, 4, 4]
# Normalize the spacing so that pixels are a distance of 1 apart
spacing = rescaled_spacing / rescaled_spacing[2]

lst = [x.split('%s_' % prefix)[1].split('_ch')[0] for x in os.listdir('%s%s' % (data_dir, sample)) if x[-4:] == '.tif']

for f in range(total_fov):
    fov = start_fov + f
    lst_temp = [int(x.split('_z')[1]) for x in lst if x.split('_z')[0] == '%s' % fov]
    total_z = max(lst_temp) + 1
    img_3d_nuclear = skio.imread("%s%s/%s_%s_z00_ch00.tif" % (data_dir, sample, prefix, fov), plugin="tifffile")[np.newaxis, ...]
    img_3d_DNAFISH = skio.imread("%s%s/%s_%s_z00_ch01.tif" % (data_dir, sample, prefix, fov), plugin="tifffile")[np.newaxis, ...]
    img_3d_laminB = skio.imread("%s%s/%s_%s_z00_ch02.tif" % (data_dir, sample, prefix, fov), plugin="tifffile")[np.newaxis, ...]
    img_3d_nuclear_seg = skio.imread("%s%s/%s_%s_z00_seg.tif" % (data_dir1, sample, prefix, fov), plugin="tifffile")[np.newaxis, ...]
    img_3d_DNAFISH_seg = skio.imread("%s%s/%s_%s_z00_ecseg1.tif" % (data_dir1, sample, prefix, fov), plugin="tifffile")[np.newaxis, ...]

    for z in range(total_z-1):
        if z < 9:
            file_name = '%s_%s_z0%s' % (prefix, fov, z+1)
        else:
            file_name = '%s_%s_z%s' % (prefix, fov, z+1)
        img_nuclear = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
        img_DNAFISH = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
        img_laminB = skio.imread("%s%s/%s_ch02.tif" % (data_dir, sample, file_name), plugin="tifffile")
        img_nuclear_seg = skio.imread("%s%s/%s_seg.tif" % (data_dir1, sample, file_name), plugin="tifffile")
        img_DNAFISH_seg = skio.imread("%s%s/%s_ecseg1.tif" % (data_dir1, sample, file_name), plugin="tifffile")
        img_3d_nuclear = np.vstack([img_3d_nuclear, img_nuclear[np.newaxis, ...]])
        img_3d_DNAFISH = np.vstack([img_3d_DNAFISH, img_DNAFISH[np.newaxis, ...]])
        img_3d_laminB = np.vstack([img_3d_laminB, img_laminB[np.newaxis, ...]])
        img_3d_nuclear_seg = np.vstack([img_3d_nuclear_seg, img_nuclear_seg[np.newaxis, ...]])
        img_3d_DNAFISH_seg = np.vstack([img_3d_DNAFISH_seg, img_DNAFISH_seg[np.newaxis, ...]])

        print(img_3d_nuclear.shape)
    tif.imwrite("%s3d/%s/%s_fov%s_DNAFISH_seg.tif" % (output_dir, sample, sample, fov), img_3d_nuclear)
    viewer = napari.Viewer()
    viewer.add_image(img_3d_nuclear, blending='additive', colormap='blue', name='nuclei', scale=original_spacing)
    viewer.add_image(img_3d_nuclear_seg, blending='additive', colormap='blue', name='nuclei seg', scale=original_spacing)
    viewer.add_image(img_3d_DNAFISH, blending='additive', colormap='green', name='DNAFISH', scale=original_spacing)
    viewer.add_image(img_3d_DNAFISH_seg, blending='additive', colormap='green', name='DNAFISH seg',
                     scale=original_spacing)
    viewer.add_image(img_3d_laminB, blending='additive', colormap='red', name='laminB', scale=original_spacing)
    # viewer.window.add_plugin_dock_widget(plugin_name='napari-animation')
    napari.run()

