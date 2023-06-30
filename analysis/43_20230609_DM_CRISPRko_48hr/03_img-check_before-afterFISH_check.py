import skimage.io as skio
import napari
import shared.image as ima
import tifffile as tif

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'E9'
topleft_target = [-30.17315988295968, -102.07718983747714]  # positive change 2nd image, negative change 1st image

img_before_GFP = skio.imread("%s/beforeFISH/%s/%s_GFP.tif" % (data_dir, sample, sample), plugin="tifffile")[:, :, 1]
img_before_mCherry = skio.imread("%s/beforeFISH/%s/%s_mCherry.tif" % (data_dir, sample, sample), plugin="tifffile")[:, :, 0]
img_after_hoechst = skio.imread("%s/afterFISH/%s/%s_hoechst.tif" % (data_dir, sample, sample), plugin="tifffile")[:, :, 2]
img_before_GFP_final, img_after_hoechst_final = ima.img_align_move(img_before_GFP, img_after_hoechst, [90, 90], topleft_target)
img_before_mCherry_final, _ = ima.img_align_move(img_before_mCherry, img_after_hoechst, [90, 90], topleft_target)

tif.imwrite("%s/beforeFISH/%s/%s_GFP_final.tif" % (output_dir, sample, sample), img_before_GFP_final)
tif.imwrite("%s/beforeFISH/%s/%s_mCherry_final.tif" % (output_dir, sample, sample), img_before_mCherry_final)
tif.imwrite("%s/afterFISH/%s/%s_hoechst_final.tif" % (output_dir, sample, sample), img_after_hoechst_final)

viewer = napari.Viewer()
viewer.add_image(img_before_GFP_final, blending='additive', colormap='green', contrast_limits=[0, 65535])
viewer.add_image(img_before_mCherry_final, blending='additive', colormap='red', contrast_limits=[0, 65535])
viewer.add_image(img_after_hoechst_final, blending='additive', colormap='blue', contrast_limits=[0, 65535])
# plt.imsave("%s%s/DM_alignment.tiff" % (output_dir, sample), dis.blending(viewer))
napari.run()