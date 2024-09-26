import nd2
import napari
import pandas as pd
import numpy as np
import shared.image as ima
import tifffile as tif

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'B3'
total_fov = 16
pixel_size = 300/2720  # uM
img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, sample))
print(img_stack.shape)
pd_loc = pd.read_csv("%s/%s/04_%s_location.txt" % (output_dir, sample, sample), na_values=['.'], sep='\t')

max_x = int((max(pd_loc['x'])-min(pd_loc['x'])+300)/pixel_size)
max_y = int((max(pd_loc['y'])-min(pd_loc['y'])+300)/pixel_size)

img_stitch = np.zeros((max_y, max_x), np.uint16)
print(img_stitch.shape)

for i in range(total_fov):
    img_hoechst = img_stack[i, :, 0, :, :]
    img_hoechst_merge = img_hoechst.max(axis=0)
    delta_x = int((pd_loc['x'][i] - min(pd_loc['x']))/pixel_size)
    delta_y = int((pd_loc['y'][i] - max(pd_loc['y']))/pixel_size)
    print(delta_x)
    print(delta_y)
    img_stitch = ima.image_paste_to(img_stitch, img_hoechst_merge, [-delta_y, delta_x])

viewer = napari.Viewer()
viewer.add_image(img_stitch, blending='additive', contrast_limits=[0, 20000])
napari.run()

tif.imwrite("%s/%s/05_%s_DNAFISH_hoechst_stitch.tif" % (output_dir, sample, sample), img_stitch)
print("DONE!")