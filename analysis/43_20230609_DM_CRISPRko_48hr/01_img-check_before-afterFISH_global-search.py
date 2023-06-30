import skimage.io as skio
import napari
import imutils
import shared.image as ima
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'E9'
interval = 50

# img
img_before_GFP = skio.imread("%s/beforeFISH/%s/%s_GFP.tif" % (data_dir, sample, sample), plugin="tifffile")[:, :, 1]
img_before_mCherry = skio.imread("%s/beforeFISH/%s/%s_mCherry.tif" % (data_dir, sample, sample), plugin="tifffile")[:, :, 0]
img_before = img_before_GFP.copy()
img_before[img_before_mCherry > img_before_GFP] = img_before_mCherry[img_before_mCherry > img_before_GFP]
img_before1 = imutils.rotate(img_before, angle=90)
img = img_before1

# img_search
img_after_hoechst = skio.imread("%s/afterFISH/%s/%s_hoechst.tif" % (data_dir, sample, sample), plugin="tifffile")[:, :, 2]
img_after_hoechst1 = imutils.rotate(img_after_hoechst, angle=90)
img_after_hoechst1_test = img_after_hoechst1[2000:5000, 2000:5000]
img_after_hoechst1_test1 = np.concatenate([np.zeros(shape=[2000, img_after_hoechst1_test.shape[1]]), img_after_hoechst1_test], axis=0)
img_after_hoechst1_test2 = np.concatenate([np.zeros(shape=[img_after_hoechst1_test1.shape[0], 2000]), img_after_hoechst1_test1], axis=1)
img_search = img_after_hoechst1_test

# search
viewer = napari.Viewer()
viewer.add_image(img_before_GFP, blending='additive', colormap='green', contrast_limits=[0, 65535])
viewer.add_image(img_before_mCherry, blending='additive', colormap='red', contrast_limits=[0, 65535])
viewer.add_image(img, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_after_hoechst1_test2, blending='additive', colormap='green', contrast_limits=[0, 65535])
shapes = viewer.add_shapes(name='Shapes', ndim=2)
napari.run()
poly_data = shapes.data[0]
topleft_target, min_ratio = ima.img_search_global(img, img_search, poly_data, interval)

# local view
img_cut = img.copy()
img_cut = img_cut[int(topleft_target[1]):int(topleft_target[1]+img_search.shape[0]), int(topleft_target[0]):int(topleft_target[0]+img_search.shape[1])]

viewer = napari.Viewer()
viewer.add_image(img_cut, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_search, blending='additive', colormap='green', contrast_limits=[0, 65535])
napari.run()

# global view
img_final, img_search_final = ima.img_align_move(img, img_search, [0, 0], topleft_target)

viewer = napari.Viewer()
viewer.add_image(img_final, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_search_final, blending='additive', colormap='green', contrast_limits=[0, 65535])
napari.run()

print([topleft_target[0]-2000, topleft_target[1]-2000])
