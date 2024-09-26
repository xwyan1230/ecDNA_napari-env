import skimage.io as skio
import napari
import cv2
import pandas as pd
import shared.image as ima
import tifffile as tif

# USAGE
# 01. square the imaging region from the screen image
# 02. add points to the topleft corner of the individual FISH images

# INPUT PARAMETERS
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"
sample = 'F8-1'
afterFISH_hoechst_pixel = 766
DNAFISH_hoechst_fov_pixel = 58.6
# 766nm for Keyence 10x
# 360nm for 512x512 at 63x
# 58.6nm for 3144x3144 at 63x (0.0765)
# 180nm for 1024x1024 at 63x (0.2349)
# 720nm for 256x256 at 63x (0.9399)

# NO NEED TO CHANGE
print(sample)
data_dir = "%sdata/" % master_folder
data_dir1 = "%sprocessed/" % master_folder
align_dir = "%salign/" % master_folder
output_dir = "%sprocessed/" % master_folder

# DO NOT CHANGE
dshape_factor = DNAFISH_hoechst_fov_pixel * 1.0 / afterFISH_hoechst_pixel
border = 100

topleft_pd = pd.DataFrame(columns=['index', 'topleft_target'])

img_after_hoechst = skio.imread("%s/%s/%s_afterFISH_hoechst_cut.tif" % (data_dir1, sample, sample), plugin="tifffile")
shape = img_after_hoechst.shape

img_screen = skio.imread("%s/screenshot/%s.png" % (data_dir, sample))
viewer = napari.Viewer()
viewer.add_image(img_screen, blending='additive', contrast_limits=[0, img_screen.max()])
shapes = viewer.add_shapes(name='Shapes', ndim=2)
napari.run()
poly_data = shapes.data[0]

img_screen_cut = img_screen[int(poly_data[0][0]):int(poly_data[2][0]), int(poly_data[0][1]):int(poly_data[1][1])]
ori_shape = img_screen_cut.shape
img_screen_cut_resize = cv2.resize(img_screen_cut, dsize=(shape[1]-2*border, shape[0]-2*border), interpolation=cv2.INTER_AREA)[:, :, 0]
_, img_screen_cut_resize = ima.img_align_move(img_after_hoechst, img_screen_cut_resize, [0, 0], [border, border])
tif.imwrite("%s/%s/%s_screen_cut.tif" % (output_dir, sample, sample), img_screen_cut_resize)

viewer = napari.Viewer()
viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_screen_cut_resize, blending='additive', contrast_limits=[0, img_screen_cut_resize.max()])
shapes = viewer.add_shapes(name='Shapes', ndim=2)
napari.run()
poly_data = shapes.data
topleft_target_lst = [x[0].tolist() for x in poly_data]

topleft_pd['topleft_target'] = topleft_target_lst
topleft_pd['index'] = topleft_pd.index
topleft_pd.to_csv('%s/%s/align_screen_topleft_locations_%s.txt' % (align_dir, sample, sample), index=False, sep='\t')

print(sample)
print("step07 find topleft from screen image Done!")
