import skimage.io as skio
import shared.segmentation as seg
import napari
import math

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220701_BRD4series/220707_BRD4_IF/"
sample = 'DM_BRD4_IF'

# cell info
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, Paul scope 2048x2048:102)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi

# IMAGING ANALYSIS
# load images
# tifs = [x for x in os.listdir("%s%s/" % (master_folder, sample)) if x[-4:] == '.tif']
# print(tifs)
img_hoechst = skio.imread("%s%s/%s_B2 Position1_RAW_ch00.tif" % (master_folder, sample, sample), plugin="tifffile")
img_IF = skio.imread("%s%s/%s_B2 Position1_RAW_ch01.tif" % (master_folder, sample, sample), plugin="tifffile")

nuclear_seg = seg.nuclear_seg(img_hoechst, local_factor=111, min_size=min_size_nuclear,
                              max_size=max_size_nuclear)

viewer = napari.view_image(img_hoechst)
viewer.add_image(nuclear_seg)
napari.run()