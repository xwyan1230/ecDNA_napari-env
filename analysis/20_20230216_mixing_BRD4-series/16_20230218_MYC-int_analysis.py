import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import math
import shared.segmentation as seg
import shared.objects as obj
import pandas as pd
import numpy as np
import os
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230216_analysis_BRD4-series/"
data_dir = "%sdata/IF/20230210_KO-series_MYC-IF/" % master_folder
output_dir = "%sfigures/IF/" % master_folder
start_fov = 1
end_fov = 5

sample = 'DM-CTCF-KO_oldMYC'
total_fov = end_fov - start_fov + 1

pixel_size = 92  # nm (sp8 confocal 3144x3144 40x oil)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
threshold = 0.8

# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi

data = pd.DataFrame(columns=['nuclear', 'FOV', 'centroid_nuclear', 'area_nuclear', 'circ_nuclear', 'mean_int_nuclear',
                             'mean_int_MYC'])

for f in range(total_fov):
    fov = start_fov + f
    file_name = '%s_Image %s_RAW' % (sample, fov)
    img_hoechst = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_MYC = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_MYC = np.concatenate([np.zeros(shape=[10, 3144]), img_MYC], axis=0)[:3144, :3144]

    # nuclear seg
    img_nuclear_seg = seg.nuclear_seg1(img_hoechst, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                       max_size=max_size_nuclear)

    img_seg = np.zeros_like(img_nuclear_seg)

    nuclear_props = regionprops(img_nuclear_seg)
    for i in range(len(nuclear_props)):
        convex_local = nuclear_props[i].convex_image
        area_ratio = nuclear_props[i].area / convex_local.sum()
        if area_ratio > threshold:
            img_seg[img_nuclear_seg == nuclear_props[i].label] = nuclear_props[i].label

    img_seg = obj.label_resort(img_seg)
    MYC_props = regionprops(img_seg, img_MYC)
    nuclear_props = regionprops(img_seg, img_hoechst)

    for i in range(len(MYC_props)):

        centroid_nuclear = MYC_props[i].centroid
        area_nuclear = MYC_props[i].area
        perimeter_nuclear = MYC_props[i].perimeter
        mean_int_nuclear = nuclear_props[i].intensity_mean
        mean_int_MYC = MYC_props[i].intensity_mean
        circ_nuclear = (4 * math.pi * area_nuclear) / (perimeter_nuclear ** 2)

        data.loc[len(data.index)] = [i, fov, centroid_nuclear, area_nuclear, circ_nuclear, mean_int_nuclear,
                                     mean_int_MYC]

    """viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
    viewer.add_image(img_seg, blending='additive', contrast_limits=[0, 1])
    napari.run()"""

    if not os.path.exists("%s%s/color_img/" % (output_dir, sample)):
        os.makedirs("%s%s/color_img/" % (output_dir, sample))

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
    plt.imsave('%s%s/color_img/%s_fov%s.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_seg, blending='additive', colormap='blue', contrast_limits=[0, 1])
    plt.imsave('%s%s/color_img/%s_fov%s_seg.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

data.to_csv('%s%s_MYC.txt' % (output_dir, sample), index=False, sep='\t')

print("DONE! ")
