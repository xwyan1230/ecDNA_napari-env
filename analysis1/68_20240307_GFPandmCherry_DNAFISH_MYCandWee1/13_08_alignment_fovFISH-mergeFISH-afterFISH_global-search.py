import skimage.io as skio
import cv2
import shared.dataframe as dat
import imutils
import pandas as pd
import shared.image as ima

# USAGE
# 01. no need to do anything

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
output_dir = "%salign/" % master_folder
name = pd.read_csv('%s/name.txt' % data_dir, na_values=['.'], sep='\t')
# treatment = name[name['sample'] == sample]['treatment'].tolist()[0]
# filename = '20240301_sp8_GFPandmCherry_MYCandWee1_4day_DNAFISH_%s_%s_%s' % (sample, treatment, sample)
filename = '20240311_GFPandmCherry_DNAFISH_F8_re-capture_%s' % sample

# DO NOT CHANGE
dshape_factor = DNAFISH_hoechst_fov_pixel * 1.0 / afterFISH_hoechst_pixel
total_fov = 16

img_after_hoechst = skio.imread("%s/%s/%s_afterFISH_hoechst_cut.tif" % (data_dir1, sample, sample), plugin="tifffile")
align_topleft = pd.read_csv('%s/%s/align_screen_topleft_locations_%s.txt' % (align_dir, sample, sample), na_values=['.'], sep='\t')
align_topleft['topleft_target'] = [dat.str_to_float(x) for x in align_topleft['topleft_target']]

topleft_search = pd.DataFrame(columns=['sample', 'fov', 'topleft_target_no', 'topleft_target', 'min_ratio'])

for fov in range(total_fov):
    print(fov+1)
    img_hoechst = skio.imread("%s/DNAFISH/%s/%s_%s_ch01.tif" % (data_dir, sample, filename, fov+1), plugin="tifffile")
    ori_shape = img_hoechst.shape
    img_hoechst_resize = cv2.resize(img_hoechst,
                                    dsize=(int(ori_shape[1] * dshape_factor), int(ori_shape[0] * dshape_factor)),
                                    interpolation=cv2.INTER_AREA)
    img_hoechst_resize1 = cv2.flip(imutils.rotate(img_hoechst_resize, angle=-90), 0)

    for i in range(len(align_topleft)):
        print("searching pos %s" % (i+1))
        topleft_target, min_ratio = ima.img_search_local(img_after_hoechst, img_hoechst_resize1, [align_topleft['topleft_target'][i][1], align_topleft['topleft_target'][i][0]], 5)
        topleft_search.loc[len(topleft_search.index)] = [sample, fov + 1, i, topleft_target, min_ratio]

topleft_search.to_csv('%s/%s/align_topleft_search_%s.txt' % (align_dir, sample, sample), index=False, sep='\t')

print(sample)
print("step08 topleft individual FOV global search Done!")