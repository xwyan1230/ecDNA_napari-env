import skimage.io as skio
import pandas as pd
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221213_analysis_DNA-FISH_mixing_exp/20221017_mixing-test_mCherry-series_after-heating/"
data_dir1 = "%sdata/" % master_folder
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'DM_mix_DM-H2B-mCherry'
file_name = '%s_RAW' % sample
img_hoechst_stack = skio.imread("%s%s_ch00.tif" % (data_dir1, file_name), plugin="tifffile")
img_mCherry_stack = skio.imread("%s%s_ch01.tif" % (data_dir1, file_name), plugin="tifffile")

data = pd.DataFrame(columns=['nuclear', 'FOV', 'mCherry_mean'])

for fov in range(img_mCherry_stack.shape[0]):
    print(fov)
    img_hoechst = img_hoechst_stack[fov, :, :]
    img_mCherry = img_mCherry_stack[fov, :, :]
    img_seg = skio.imread("%spdf/01_data_tif/%s/seg_tif/%s_%s_seg.tif" % (data_dir2, sample, sample, fov), plugin="tifffile")

    # measure
    nuclear_props = regionprops(label(img_seg), img_mCherry)
    for j in range(len(nuclear_props)):
        data.loc[len(data.index)] = [nuclear_props[j].label, fov, nuclear_props[j].intensity_mean]
        print(nuclear_props[j].area)
        print(nuclear_props[j].intensity_mean)

# data.to_csv('%s%s.txt' % (output_dir, sample), index=False, sep='\t')
print("DONE!")
