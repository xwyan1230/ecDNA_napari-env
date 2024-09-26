import skimage.io as skio
import napari
import numpy as np
import pandas as pd
import shared.image as ima
import nd2
import shared.dataframe as dat
from skimage.measure import label, regionprops
from skimage.morphology import medial_axis

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'G2'
total_fov = 16
dshape_factor = 0.145
pixel_size = 300 / 2720  # uM
local_size = 200
img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, sample))

data = pd.DataFrame(columns=['sample', 'fov', 'label_nuclear',
                             'DNAFISH_seg_label', 'int_r_to_edge', 'int_relative_r'])

pd_seg_z = pd.read_csv('%s/%s/10_%s_seg_z.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask.flatten()) if e != 0]
    out = list(map(img.flatten().__getitem__, index))
    return out


for fov in range(total_fov):
    print("%s/%s" % (fov + 1, total_fov))
    img_hoechst = img_stack[fov, :, 0, :, :]
    img_DNAFISH = img_stack[fov, :, 1, :, :]
    seg_z = pd_seg_z['seg_z'][fov]
    img_hoechst_seg_z = img_hoechst[seg_z]
    img_DNAFISH_seg_z = img_DNAFISH[seg_z]

    img_seg_z = skio.imread("%s/%s/13_seg_new_tif/%s_%s_seg_z_new.tif" % (output_dir, sample, sample, fov),
                            plugin="tifffile")
    img_seg_ecDNA = skio.imread("%s/%s/17_seg_DNAFISH_tif/%s_%s_seg_DNAFISH.tif" % (output_dir, sample, sample, fov),
                            plugin="tifffile")

    nuclear_props = regionprops(label(img_seg_z), img_seg_z)

    for i in range(len(nuclear_props)):
        print("Analyzing %s, fov %s, nuclear %s/%s" % (sample, fov, i + 1, len(nuclear_props)))
        original_centroid_nuclear = nuclear_props[i].centroid
        label_nuclear = nuclear_props[i].intensity_mean
        position = ima.img_local_position(img_seg_z, original_centroid_nuclear, local_size)
        local_nuclear_seg = ima.img_local_seg(img_seg_z, position, nuclear_props[i].intensity_mean)
        local_nuclear = img_hoechst_seg_z.copy()[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH_seg_z.copy()[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg = img_seg_ecDNA.copy()[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg[local_nuclear_seg == 0] = 0
        local_nuclear_props = regionprops(label(local_nuclear_seg), local_nuclear)

        # radial measurements
        local_nuclear_centroid = local_nuclear_props[0].centroid
        _, local_edge_distance_map = medial_axis(local_nuclear_seg, return_distance=True)
        local_centroid_distance_map = ima.distance_map_from_point(local_nuclear_seg, local_nuclear_centroid)
        local_centroid_distance_map[local_nuclear_seg == 0] = 0
        local_edge_distance_map[local_nuclear_seg == 0] = -1
        local_relative_r_map = local_centroid_distance_map / (
                local_centroid_distance_map + local_edge_distance_map)
        DNAFISH_seg_label = img_to_pixel_int(local_nuclear_seg, local_DNAFISH_seg)
        int_r_to_edge = img_to_pixel_int(local_nuclear_seg, local_edge_distance_map)
        int_relative_r = img_to_pixel_int(local_nuclear_seg, local_relative_r_map)

        data.loc[len(data.index)] = [sample, fov, label_nuclear, DNAFISH_seg_label, int_r_to_edge, int_relative_r]

data.to_csv('%s/%s/18_%s_radial.txt' % (output_dir, sample, sample), index=False, sep='\t')
print("DONE!")
