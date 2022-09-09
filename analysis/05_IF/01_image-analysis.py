import skimage.io as skio
import shared.segmentation as seg
from skimage.measure import regionprops, label
import napari
import pandas as pd
import math

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220922_EdU_test/"
sample = 'CNTRL'
prefix = '220831_'
postfix = '_1_RAW'
sample_name = '%s%s%s' % (prefix, sample, postfix)

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
img_hoechst_stack = skio.imread("%s%s_ch00.tif" % (master_folder, sample_name), plugin="tifffile")
img_IF_stack = skio.imread("%s%s_ch01.tif" % (master_folder, sample_name), plugin="tifffile")

data = pd.DataFrame(columns=['sample', 'FOV', 'nuclear', 'mean_int_IF'])

total_fov = img_hoechst_stack.shape[0]

for fov in range(total_fov):
    print("Analyzing %s, %s/%s..." % (sample, fov+1, total_fov))
    nuclear_seg = seg.nuclear_seg(img_hoechst_stack[fov], local_factor=111, min_size=min_size_nuclear,
                                  max_size=max_size_nuclear)
    nuclear_props = regionprops(label(nuclear_seg), img_IF_stack[fov])
    for n in range(len(nuclear_props)):
        data.loc[len(data.index)] = [sample, fov, n, nuclear_props[n].intensity_mean]

    """viewer = napari.view_image(img_hoechst_stack[fov])
    viewer.add_image(nuclear_seg)
    napari.run()"""

data.to_csv('%s%s.txt' % (master_folder, sample), index=False, sep='\t')

print("DONE!")