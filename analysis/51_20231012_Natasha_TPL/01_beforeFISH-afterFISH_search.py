import skimage.io as skio
import napari
import cv2
import imutils
import tifffile as tif
import shared.image as ima

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230724_analysis_Natasha_colcemid/TPL/mCh-DMSO_GFP-TPL/"
output_dir = "%s05_figures/" % master_folder

sample = 'mCh-DMSO_GFP-TPL'
interval = 50
dshape_factor = 0.471  # 766nm for Keyence 10x, 360nm for 512x512 at 63x, 58.6 for 3144x3144 at 63x (0.0765), 180 for 1024x1024 at 63x (0.2349)
frame = 4

# img
img_before = skio.imread("%s/01_beforeFISH/DAPI.tif" % master_folder, plugin="tifffile")[:, :, 2]
img_before_red = skio.imread("%s/01_beforeFISH/mCherry.tif" % master_folder, plugin="tifffile")[:, :, 0]
img_before_green = skio.imread("%s/01_beforeFISH/GFP.tif" % master_folder, plugin="tifffile")[:, :, 1]
img = imutils.rotate(img_before, angle=-1)
img_before_red_rotate = imutils.rotate(img_before_red, angle=-1)
img_before_green_rotate = imutils.rotate(img_before_green, angle=-1)
img_after = skio.imread("%s/02_afterFISH/230606_mChDMSO_GFPTPL_mChDMSO_GFPTPL_overview_4_Merged_ch00.tif" % master_folder, plugin="tifffile")
ori_shape = img_after.shape
img_search = cv2.resize(img_after, dsize=(int(ori_shape[1]*dshape_factor), int(ori_shape[0]*dshape_factor)), interpolation=cv2.INTER_AREA)

"""# global search
viewer = napari.Viewer()
viewer.add_image(img, blending='additive', colormap='blue', contrast_limits=[0, img.max()])
viewer.add_image(img_search, blending='additive', contrast_limits=[0, img_search.max()])
shapes = viewer.add_shapes(name='Shapes', ndim=2)
napari.run()
poly_data = shapes.data[0]
topleft_target, min_ratio = ima.img_search_global(img, img_search, poly_data, interval)"""

topleft_target = [2932, 2670]

# local view
print(topleft_target)
print(sample)
img_cut = img.copy()
img_cut = img_cut[int(topleft_target[1]):int(topleft_target[1]+img_search.shape[0]), int(topleft_target[0]):int(topleft_target[0]+img_search.shape[1])]
tif.imwrite("%s/05_figures/%s_img_hoechst_cut.tif" % (master_folder, frame), img_cut)
img_red_cut = img_before_red_rotate.copy()[int(topleft_target[1]):int(topleft_target[1]+img_search.shape[0]), int(topleft_target[0]):int(topleft_target[0]+img_search.shape[1])]
img_green_cut = img_before_green_rotate.copy()[int(topleft_target[1]):int(topleft_target[1]+img_search.shape[0]), int(topleft_target[0]):int(topleft_target[0]+img_search.shape[1])]
tif.imwrite("%s/05_figures/%s_img_red_cut.tif" % (master_folder, frame), img_red_cut)
tif.imwrite("%s/05_figures/%s_img_green_cut.tif" % (master_folder, frame), img_green_cut)

viewer = napari.Viewer()
viewer.add_image(img_cut, blending='additive', colormap='blue', contrast_limits=[0, img_cut.max()])
viewer.add_image(img_search, blending='additive', colormap='green', contrast_limits=[0, img_search.max()])
napari.run()

# global view
img_final, img_search_final = ima.img_align_move(img, img_search, [0, 0], topleft_target)

viewer = napari.Viewer()
viewer.add_image(img_final, blending='additive', colormap='blue', contrast_limits=[0, img_final.max()])
viewer.add_image(img_search_final, blending='additive', colormap='green', contrast_limits=[0, img_search_final.max()])
napari.run()