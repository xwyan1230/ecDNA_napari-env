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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230711_analysis_H3K4me3_NPC_PC3/"
data_dir = "%sdata/raw/" % master_folder
output_dir = "%sdata/" % master_folder

# sample = 'Colo320DM'
# sample = 'NPC'
# sample = 'TPR'
# sample = 'Nup153'
sample = 'PC3DM'
# prefix = '20230616_H3K4me3_DM_H3K4me3'
# prefix = '20230608_DM_IF-NPC_acidFISH-MYC_ColoDM_IF-NPC_acidFISH-MYC_100x-2048'
# prefix = '20230616_DM_TPR-Nup98_DM_TPR'
# prefix = '20230608_DM_IF-NPC_acidFISH-MYC_ColoDM_IF-Nup153_acidFISH-MYC_100x-2048_z-point3um'
prefix = '20230616_PC3_heat_lamin_PC3DM_heat_lamin'
total_fov = 3
start_fov = 1

# The microscope reports the following spacing (in µm)
original_spacing = np.array([0.3, 0.05679, 0.05679])
# We downsampled each slice 4x to make the data smaller
rescaled_spacing = original_spacing * [1, 4, 4]
# Normalize the spacing so that pixels are a distance of 1 apart
spacing = rescaled_spacing / rescaled_spacing[2]

lst = [x.split('%s_' % prefix)[1].split('_ch')[0] for x in os.listdir('%s%s' % (data_dir, sample)) if x[-4:] == '.tif']
print(lst)

for f in range(total_fov):
    fov = start_fov + f
    lst_temp = [int(x.split('_z')[1]) for x in lst if x.split('_z')[0] == '%s' % fov]
    total_z = max(lst_temp) + 1
    img_3d_nuclear = skio.imread("%s%s/%s_%s_z00_ch00.tif" % (data_dir, sample, prefix, fov), plugin="tifffile")[np.newaxis, ...]
    img_3d_DNAFISH = skio.imread("%s%s/%s_%s_z00_ch01.tif" % (data_dir, sample, prefix, fov), plugin="tifffile")[np.newaxis, ...]
    img_3d_laminB = skio.imread("%s%s/%s_%s_z00_ch02.tif" % (data_dir, sample, prefix, fov), plugin="tifffile")[np.newaxis, ...]
    for z in range(total_z-1):
        if z < 9:
            file_name = '%s_%s_z0%s' % (prefix, fov, z+1)
        else:
            file_name = '%s_%s_z%s' % (prefix, fov, z+1)
        img_nuclear = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
        img_DNAFISH = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
        img_laminB = skio.imread("%s%s/%s_ch02.tif" % (data_dir, sample, file_name), plugin="tifffile")
        img_3d_nuclear = np.vstack([img_3d_nuclear, img_nuclear[np.newaxis, ...]])
        img_3d_DNAFISH = np.vstack([img_3d_DNAFISH, img_DNAFISH[np.newaxis, ...]])
        img_3d_laminB = np.vstack([img_3d_laminB, img_laminB[np.newaxis, ...]])
        print(img_3d_nuclear.shape)
    # img_4d = np.vstack([img_3d_nuclear[np.newaxis, ...], img_3d_DNAFISH[np.newaxis, ...], img_3d_laminB[np.newaxis, ...]])
    # print(img_4d.shape)
    """if not os.path.exists("%s3d/%s/" % (output_dir, sample)):
        os.makedirs("%s3d/%s/" % (output_dir, sample))
    tif.imwrite("%s3d/%s/%s_fov%s_nuclear.tif" % (output_dir, sample, sample, fov), img_3d_nuclear)
    tif.imwrite("%s3d/%s/%s_fov%s_DNAFISH.tif" % (output_dir, sample, sample, fov), img_3d_DNAFISH)
    tif.imwrite("%s3d/%s/%s_fov%s_laminB.tif" % (output_dir, sample, sample, fov), img_3d_laminB)"""

    viewer = napari.Viewer()
    viewer.add_image(img_3d_nuclear, blending='additive', colormap='blue', name='nuclei', scale=original_spacing)
    viewer.add_image(img_3d_DNAFISH, blending='additive', colormap='green', name='DNAFISH', scale=original_spacing)
    viewer.add_image(img_3d_laminB, blending='additive', colormap='red', name='laminB', scale=original_spacing)
    # viewer.window.add_plugin_dock_widget(plugin_name='napari-animation')
    napari.run()

