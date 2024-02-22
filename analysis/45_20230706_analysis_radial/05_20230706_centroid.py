import skimage.io as skio
import napari
import imutils
import shared.image as ima
from skimage.measure import label, regionprops
import numpy as np
import pandas as pd
# from napari_animation import Animation
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230706_analysis_radial/"
data_dir = "%sdata/raw/" % master_folder
data_dir1 = "%sdata/seg_tif/" % master_folder
output_dir = "%sfigures/" % master_folder

# sample = 'Colo320DM_acidFISH_lamin_3d'
sample = 'Colo320HSR_acidFISH_lamin_3d'
# prefix = '20230614_acidFISH_lamin_ColoDM_DM'
prefix = '20230614_acidFISH_lamin_ColoHSR_HSR'
total_fov = 6
start_fov = 1

# The microscope reports the following spacing (in µm)
original_spacing = np.array([0.3, 0.05679, 0.05679])
# We downsampled each slice 4x to make the data smaller
rescaled_spacing = original_spacing * [1, 4, 4]
# Normalize the spacing so that pixels are a distance of 1 apart
spacing = rescaled_spacing / rescaled_spacing[2]

lst = [x.split('%s_' % prefix)[1].split('_ch')[0] for x in os.listdir('%s%s' % (data_dir, sample)) if x[-4:] == '.tif']

data = pd.DataFrame(columns=['FOV', 'z', 'label', 'area', 'centroid_x', 'centroid_y'])

for f in range(total_fov):
    fov = start_fov + f
    lst_temp = [int(x.split('_z')[1]) for x in lst if x.split('_z')[0] == '%s' % fov]
    total_z = max(lst_temp) + 1
    for z in range(total_z):
        if z < 10:
            file_name = '%s_%s_z0%s' % (prefix, fov, z)
        else:
            file_name = '%s_%s_z%s' % (prefix, fov, z)
        img_nuclear_seg = skio.imread("%s%s/%s_seg.tif" % (data_dir1, sample, file_name), plugin="tifffile")
        nuclear_props = regionprops(img_nuclear_seg)
        for i in range(len(nuclear_props)):
            data.loc[len(data.index)] = [fov, z, nuclear_props[i].label, nuclear_props[i].area, nuclear_props[i].centroid[0],
                                         nuclear_props[i].centroid[1]]

if not os.path.exists("%s%s/" % (output_dir, sample)):
    os.makedirs("%s%s/" % (output_dir, sample))
data.to_csv('%s%s/%s_centroid.txt' % (output_dir, sample, sample), index=False, sep='\t')
print("DONE!")



