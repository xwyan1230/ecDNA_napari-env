import skimage.io as skio
import napari
import tifffile as tif
import shared.display as dis
import matplotlib.pyplot as plt

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/papers/2024_Wee1/Figures/data/17_FigS2c_cytokinesis/"

sample = 'BDP5290_5uM_48hr_2_XY40_1'
img_hoechst = skio.imread("%s/%s.tif" % (master_folder, sample), plugin="tifffile")[:, :, 2]
print(img_hoechst.shape)
width = 520  # 200um

viewer = napari.Viewer()
viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 40000])
shapes = viewer.add_shapes(name='Shapes', ndim=2)
napari.run()
center_data_x = int(shapes.data[0][0][0])
center_data_y = int(shapes.data[0][0][1])

img_crop = img_hoechst.copy()
img_crop = img_crop[int(center_data_x-width/2):int(center_data_x+width/2), int(center_data_y-width/2):int(center_data_y+width/2)]

tif.imwrite("%s/%s_crop.tif" % (master_folder, sample), img_crop)
img_crop = skio.imread("%s/%s_crop.tif" % (master_folder, sample), plugin="tifffile")
print(img_crop.shape)

viewer = napari.Viewer()
viewer.add_image(img_crop, blending='additive', colormap='blue', contrast_limits=[0, 40000])
plt.imsave("%s/%s_crop.tiff" % (master_folder, sample), dis.blending(viewer))
napari.run()


