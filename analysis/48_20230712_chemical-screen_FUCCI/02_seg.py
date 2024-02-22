import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import tifffile as tif
from skimage.measure import regionprops
import numpy as np
import shared.segmentation as seg
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230712_analysis_chemical-screen_FUCCI/"
data_dir = "%sdata/raw/" % master_folder
output_dir = "%sdata/" % master_folder
output_dir1 = "%sfigures/" % master_folder

plate = 'HSR_FUCCI_2hr'
sample_lst = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'C11', 'C10', 'C9', 'C8', 'C7', 'C6',
              'C5', 'C4', 'C3', 'C2', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'E11', 'E10',
              'E9', 'E8', 'E7', 'E6', 'E5', 'E4', 'E3', 'E2', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10',
              'F11', 'G11', 'G10', 'G9', 'G8', 'G7', 'G6', 'G5', 'G4', 'G3', 'G2']

# set parameters
otsu_factor = 1.2
circ_threshold = 0.8
min_size = 400
max_size = 2000

start_sample = 'F10'
start_index = sample_lst.index(start_sample)

for s in range(len(sample_lst)):
    sample = sample_lst[s+start_index]
    data = pd.DataFrame(columns=['plate', 'well', 'nuclear', 'FOV', 'centroid_nuclear', 'area_nuclear', 'mean_int_nuclear', 'mean_int_green', 'mean_int_red'])
    print(sample)
    pos = sample_lst.index(sample)
    print(pos)
    for i in range(4):
        df = pd.DataFrame()
        file_num = 4*pos+i+1
        if file_num < 10:
            file_folder = 'XY0%s' % file_num
        else:
            file_folder = 'XY%s' % file_num
        img_hoechst = skio.imread("%s%s/%s/Image_%s_CH1.tif" % (data_dir, plate, file_folder, file_folder), plugin="tifffile")[..., 2]  # 1536, 2048
        img_green = skio.imread("%s%s/%s/Image_%s_CH2.tif" % (data_dir, plate, file_folder, file_folder), plugin="tifffile")[..., 1]
        img_red = skio.imread("%s%s/%s/Image_%s_CH3.tif" % (data_dir, plate, file_folder, file_folder), plugin="tifffile")[..., 0]

        # nuclear seg
        img_nuclear_seg = seg.cell_seg_fluorescent(img_hoechst, otsu_factor=otsu_factor, max_size=max_size,
                                                   min_size=min_size, circ_thresh=circ_threshold).astype(int)

        nuclear_props = regionprops(img_nuclear_seg, img_hoechst)
        green_props = regionprops(img_nuclear_seg, img_green)
        red_props = regionprops(img_nuclear_seg, img_red)

        df['plate'] = [plate] * len(nuclear_props)
        df['well'] = [sample] * len(nuclear_props)
        df['nuclear'] = [nuclear_props[x].label for x in range(len(nuclear_props))]
        df['FOV'] = [i] * len(nuclear_props)
        df['centroid_nuclear'] = [nuclear_props[x].centroid for x in range(len(nuclear_props))]
        df['area_nuclear'] = [nuclear_props[x].area for x in range(len(nuclear_props))]
        df['mean_int_nuclear'] = [nuclear_props[x].intensity_mean for x in range(len(nuclear_props))]
        df['mean_int_green'] = [green_props[x].intensity_mean for x in range(len(nuclear_props))]
        df['mean_int_red'] = [red_props[x].intensity_mean for x in range(len(nuclear_props))]
        data = pd.concat([data, df], axis=0)

        if not os.path.exists("%s/seg_tif/%s/" % (output_dir, plate)):
            os.makedirs("%s/seg_tif/%s/" % (output_dir, plate))
        tif.imwrite("%s/seg_tif/%s/%s_%s_seg.tif" % (output_dir, plate, sample, file_folder), img_nuclear_seg)

        if not os.path.exists("%s/color_img/%s/" % (output_dir, plate)):
            os.makedirs("%s/color_img/%s/" % (output_dir, plate))

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 255])
        viewer.add_image(img_nuclear_seg, blending='additive')
        plt.imsave("%s/color_img/%s/%s_%s_seg.tiff" % (output_dir, plate, sample, file_folder), dis.blending(viewer))
        viewer.close()
    if not os.path.exists("%s%s/txt/" % (output_dir1, plate)):
        os.makedirs("%s%s/txt/" % (output_dir1, plate))
    data.to_csv('%s%s/txt/%s.txt' % (output_dir1, plate, sample), index=False, sep='\t')
print("DONE!")

