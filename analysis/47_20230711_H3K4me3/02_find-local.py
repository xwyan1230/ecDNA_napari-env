import skimage.io as skio
import napari
import imutils
import shared.image as ima
import numpy as np
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
# sample = 'Nup153'
# prefix = '20230616_H3K4me3_DM_H3K4me3'
# prefix = '20230608_DM_IF-NPC_acidFISH-MYC_ColoDM_IF-NPC_acidFISH-MYC_100x-2048'
prefix = '20230616_DM_TPR-Nup98_DM_TPR'
# prefix = '20230608_DM_IF-NPC_acidFISH-MYC_ColoDM_IF-Nup153_acidFISH-MYC_100x-2048_z-point3um'
total_fov = 2
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
    fov = 1
    lst_temp = [int(x.split('_z')[1]) for x in lst if x.split('_z')[0] == '%s' % fov]
    total_z = max(lst_temp) + 1
    img_3d_nuclear = skio.imread("%s%s/%s_%s_z17_ch00.tif" % (data_dir, sample, prefix, fov), plugin="tifffile")[np.newaxis, ...]
    img_3d_DNAFISH = skio.imread("%s%s/%s_%s_z17_ch01.tif" % (data_dir, sample, prefix, fov), plugin="tifffile")[np.newaxis, ...]
    img_3d_laminB = skio.imread("%s%s/%s_%s_z17_ch02.tif" % (data_dir, sample, prefix, fov), plugin="tifffile")[np.newaxis, ...]
    viewer = napari.Viewer()
    viewer.add_image(img_3d_nuclear, blending='additive', colormap='blue', name='nuclei')
    viewer.add_image(img_3d_DNAFISH, blending='additive', colormap='green', name='DNAFISH')
    viewer.add_image(img_3d_laminB, blending='additive', colormap='red', name='laminB')
    shapes = viewer.add_shapes(name='Shapes', ndim=2)
    napari.run()
    shapes_layer = viewer.layers['Shapes']
    print(shapes.data)

    # viewer.window.add_plugin_dock_widget(plugin_name='napari-animation')

