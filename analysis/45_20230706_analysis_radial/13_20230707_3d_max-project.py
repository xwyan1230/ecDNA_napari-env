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
output_dir = "%sdata/" % master_folder

# sample = 'Colo320DM_acidFISH_lamin_3d'
sample = 'Colo320HSR_acidFISH_lamin_3d'
# prefix = '20230614_acidFISH_lamin_ColoDM_DM'
prefix = '20230614_acidFISH_lamin_ColoHSR_HSR'
total_fov = 6
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
    img_3d_DNAFISH_proj = np.max(img_3d_DNAFISH, axis=0)
    if not os.path.exists("%smax_proj/%s/" % (output_dir, sample)):
        os.makedirs("%smax_proj/%s/" % (output_dir, sample))
    tif.imwrite("%smax_proj/%s/fov%s_DNAFISH_maxproj.tif" % (output_dir, sample, fov), img_3d_DNAFISH_proj)
    viewer = napari.Viewer()
    viewer.add_image(img_3d_DNAFISH_proj, blending='additive', colormap='green', name='DNAFISH')
    # viewer.window.add_plugin_dock_widget(plugin_name='napari-animation')
    napari.run()

