import skimage.io as skio
import napari
import imutils
import shared.image as ima
import shared.display as dis
import matplotlib.pyplot as plt
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/" % master_folder
output_dir = "%salignment/" % master_folder

sample = 'E9'
topleft_target = [1964.8268401170403, 1912.9228101625229]
interval = 5

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
img_search = img_after_hoechst1_test

# local_search
topleft_target, min_ratio = ima.img_search_local(img, img_search, topleft_target, interval)

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
if not os.path.exists("%s%s/" % (output_dir, sample)):
    os.makedirs("%s%s/" % (output_dir, sample))
plt.imsave("%s%s/%s_before-afterFISH.tiff" % (output_dir, sample, sample), dis.blending(viewer))
napari.run()

print([topleft_target[0]-2000, topleft_target[1]-2000])
