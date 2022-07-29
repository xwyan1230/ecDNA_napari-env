import pandas as pd
import skimage.io as skio
import shared.dataframe as dat
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220526_flowFISH_topHits_screen/"
sample = 'C4'
raw_folder = '01_raw'
seg_folder = '02_seg'
# cell info
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, Paul scope 2048x2048:102)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
z_size = 500  # nm (Paul scope)
local_size = 100
rmax = 100

# LOAD Z FILE
data_z = pd.read_csv('%s%s/%s/%s_z.txt' % (master_folder, sample[0], sample[1:], sample), na_values=['.'], sep='\t')
data_z['centroid_nuclear'] = [dat.str_to_float(data_z['centroid_nuclear'][i]) for i in range(len(data_z))]

data_z_fov0 = data_z[data_z['FOV'] == 0].copy()
z_max = max(data_z_fov0['z'].tolist())
z_min = min(data_z_fov0['z'].tolist())
data_z_fov0['modified_z'] = (data_z['z'] - z_min) * 100 / (z_max-z_min)

im_z_stack_nuclear = skio.imread("%s%s/%s/%s/R%s_RAW_ch00.tif" % (master_folder, sample[0], sample[1:], raw_folder, 1),
                                 plugin="tifffile")

viewer = napari.view_image(im_z_stack_nuclear)
points = data_z_fov0['centroid_nuclear'].tolist()
point_properties = {
    'label': data_z_fov0['label_nuclear'].tolist(),
    'z': data_z_fov0['z'].tolist(),
    'modified_z': data_z_fov0['modified_z'].tolist(),
}

points_layer = viewer.add_points(
    points,
    properties=point_properties,
    face_color='modified_z',
    face_colormap='PiYG',
    edge_color='modified_z',
    edge_colormap='PiYG',
    size=30,
)
napari.run()

