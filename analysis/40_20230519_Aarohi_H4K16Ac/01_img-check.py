import skimage.io as skio
import napari
import pandas as pd

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230519_analysis_Aarohi_H4K16Ac-test/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'H4K16Ac_untreated'
img_hoechst = skio.imread("%s/%s/%s_Hoechst.tif" % (data_dir, sample, sample), plugin="tifffile")[:, :, 2]

data = pd.DataFrame(columns=['cell', 'center_x', 'center_y'])

viewer = napari.Viewer()
viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
shapes = viewer.add_shapes(name='Shapes', ndim=2)
napari.run()
shapes_layer = viewer.layers['Shapes']

for i in range(len(shapes.data)):
    center = shapes.data[i][0]
    data.loc[len(data.index)] = [int(i), center[0], center[1]]

data.to_csv('%s%s_mitosis_location.txt' % (output_dir, sample), index=False, sep='\t')
print("DONE!")
