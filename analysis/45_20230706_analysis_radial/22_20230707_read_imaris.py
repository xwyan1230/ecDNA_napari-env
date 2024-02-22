from imaris_ims_file_reader.ims import ims
import numpy as np
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230706_analysis_radial/"
data_dir = "%sdata/imaris/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'Colo320DM_acidFISH_lamin_3d'
# sample = 'Colo320HSR_acidFISH_lamin_3d'
original_spacing = np.array([0.3, 0.05679, 0.05679])

test = ims('%s20230614_acidFISH_lamin_ColoDM_[ii0_DM_1_Image_1]_1.ims' % data_dir)[0]
print(test.shape)

img_3d_nuclear = test[0, ...]
img_3d_DNAFISH = test[1, ...]
img_3d_laminB = test[2, ...]
img_3d_nuclear_surface = test[3, ...]

viewer = napari.Viewer()
viewer.add_image(img_3d_nuclear, blending='additive', colormap='blue', name='nuclei', scale=original_spacing)
viewer.add_image(img_3d_DNAFISH, blending='additive', colormap='green', name='DNAFISH', scale=original_spacing)
viewer.add_image(img_3d_laminB, blending='additive', colormap='red', name='laminB', scale=original_spacing)
viewer.add_image(img_3d_nuclear_surface, blending='additive', name='surface', scale=original_spacing)
# viewer.window.add_plugin_dock_widget(plugin_name='napari-animation')
napari.run()