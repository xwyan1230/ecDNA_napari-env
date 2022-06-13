import pandas as pd
import shared.image as img
import skimage.io as skio
import shared.dataframe as dat
import tifffile as tif
import matplotlib.pyplot as plt
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220325_Natasha_THZ1/"
sample = '15hr_JQ1'
raw_folder = '01_raw'
seg_folder = '02_seg'
save_folder = '04_DNAFISH_scale'
total_fov = 31
total_z = 19
# cell info
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, Paul scope 2048x2048:102)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
z_size = 500  # nm (Paul scope)
local_size = 100

# LOAD Z FILE
data_z = pd.read_csv('%s%s/%s_z.txt' % (master_folder, sample, sample), na_values=['.'], sep='\t')
data_z['centroid_nuclear'] = [dat.str_to_float(data_z['centroid_nuclear'][i]) for i in range(len(data_z))]

# IMAGING ANALYSIS
for i in range(len(data_z)):
    print("Analyzing nucleus %s/%s" % (i+1, len(data_z)))
    fov = data_z['FOV'][i]
    z_current = data_z['z'][i]
    label_nuclear = data_z['label_nuclear'][i]
    original_centroid_nuclear = data_z['centroid_nuclear'][i]
    # load images
    im_z_stack_DNAFISH = skio.imread("%s%s/%s/%s_RAW_ch01_fov%s.tif" % (master_folder, sample, raw_folder, sample, fov),
                                     plugin="tifffile")
    im_z_stack_seg_convex = skio.imread("%s%s/%s/%s_seg_fov%s.tif" % (master_folder, sample, seg_folder, sample, fov),
                                        plugin="tifffile")
    # get images for given z
    img_nuclear_seg_convex = im_z_stack_seg_convex[z_current]
    img_DNAFISH = im_z_stack_DNAFISH[z_current]
    # get local images
    position = img.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear, local_size)
    local_DNAFISH = img_DNAFISH.copy()
    local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
    plt.imsave('%s%s/%s/%s_DNAFISH_fov%s_z%s_i%s.tiff' % (master_folder, sample, save_folder, sample, fov, z_current,
                                                         label_nuclear), local_DNAFISH, vmin=0, vmax=5000)
    # tif.imwrite('%s%s/%s/%s_DNAFISH_fov%s_z%s_i%s.tif' % (master_folder, sample, save_folder, sample, fov, z_current,
    #                                                       label_nuclear), local_DNAFISH)

print("DONE!")
