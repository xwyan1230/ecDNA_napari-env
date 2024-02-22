import skimage.io as skio
import napari
import cv2
import shared.image as ima
import imutils
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

row = 'G'
sample = 'G9'
dshape_factor = 0.2349  # 766nm for Keyence 10x, 360nm for 512x512 at 63x, 58.6 for 3144x3144 at 63x (0.0765), 180 for 1024x1024 at 63x (0.2349)
interval = 50

# img
# filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_scan_1024_Merging_001_ch00' % (row, sample)
filename = '20230602_DM_CRISPRko_48hr_DNAFISH_%s_%s_scan_1024_Merging_001_ch00' % (row, sample)
img_after_hoechst = skio.imread("%s/afterFISH/%s/%s_hoechst_final.tif" % (data_dir, sample, sample), plugin="tifffile")
img = img_after_hoechst

# img_search
img_hoechst_DNAFISH = skio.imread("%s/FISH/%s/%s.tif" % (data_dir, sample, filename), plugin="tifffile")
ori_shape = img_hoechst_DNAFISH.shape
img_hoechst_DNAFISH1 = cv2.resize(img_hoechst_DNAFISH, dsize=(int(ori_shape[1]*dshape_factor), int(ori_shape[0]*dshape_factor)), interpolation=cv2.INTER_AREA)
img_search = img_hoechst_DNAFISH1

# referene images
img_ref_screen = skio.imread("%s/alignment/E5/E5_screen.png" % master_folder)
img_ref_screen1 = cv2.resize(img_ref_screen, dsize=(int(img_ref_screen.shape[1]*55), int(img_ref_screen.shape[0]*55)), interpolation=cv2.INTER_AREA)
img_screen = skio.imread("%s/alignment/%s/%s_screen.png" % (master_folder, sample, sample))
img_screen1 = cv2.resize(img_screen, dsize=(int(img_screen.shape[1]*55), int(img_screen.shape[0]*55)), interpolation=cv2.INTER_AREA)
img_ref = skio.imread("%s/alignment/E5/E5_afterFISH-FISH.tiff" % master_folder)

# global search
viewer = napari.Viewer()
viewer.add_image(img, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_screen1, blending='additive', contrast_limits=[0, img_screen.max()])
viewer.add_image(img_ref, blending='additive', contrast_limits=[0, img_ref.max()])
viewer.add_image(img_ref_screen1, blending='additive', contrast_limits=[0, img_ref_screen.max()])
shapes = viewer.add_shapes(name='Shapes', ndim=2)
napari.run()
poly_data = shapes.data[0]
topleft_target, min_ratio = ima.img_search_global(img, img_search, poly_data, interval)

# local view
print(topleft_target)
print(sample)
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

