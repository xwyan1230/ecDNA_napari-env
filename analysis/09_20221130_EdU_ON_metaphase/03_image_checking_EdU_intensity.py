import skimage.io as skio
import napari
import pandas as pd
import shared.dataframe as dat
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221208_analysis_20221130_EdU_ON_metaphase/figures/"
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221208_analysis_20221130_EdU_ON_metaphase/figures/"
sample = 'point5uM_EdU'

# load and format data
df = pd.read_csv(("%s%s.txt" % (master_folder, sample)), na_values=['.'], sep='\t')
feature = ['EdU_ind_mean_int', 'EdU_ind_area', 'MYC_ind_mean_int']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]
df['EdU_centroid'] = [dat.str_to_list_of_float(df['EdU_centroid'][i], 2) for i in range(len(df))]

# filter data
df = df[df['n'] > 10].copy().reset_index(drop=True)

# data analysis
for i in range(len(df)):
    fov = df['FOV'][i]
    nuclear = df['nuclear'][i]
    points = df['EdU_centroid'][i]
    print(nuclear)

    cell_hoechst = skio.imread("%sindividual_cell_tif/%s_%s_%s_hoechst.tif" % (data_dir, sample, fov, nuclear), plugin="tifffile")
    cell_FISH = skio.imread("%sindividual_cell_tif/%s_%s_%s_FISH.tif" % (data_dir, sample, fov, nuclear), plugin="tifffile")
    cell_EdU = skio.imread("%sindividual_cell_tif/%s_%s_%s_EdU.tif" % (data_dir, sample, fov, nuclear), plugin="tifffile")

    viewer = napari.Viewer()
    viewer.add_image(cell_FISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(cell_EdU, blending='additive', colormap='magenta', contrast_limits=[0, 65535])

    features = {'EdU_mean': df['EdU_ind_mean_int'][i],
                'EdU_normalized_mean': list(np.array(df['EdU_ind_mean_int'][i])/np.array(df['MYC_ind_mean_int'][i]))}

    text = {'string': '{EdU_normalized_mean:.2f}', 'size': 10, 'color': 'white', 'translation': np.array([-5, 5])}
    points_layer = viewer.add_points(points, features=features, text=text, size=5, edge_width=1,
                                     edge_width_is_relative=False, edge_color='EdU_normalized_mean', edge_colormap='PiYG',
                                     face_color='EdU_normalized_mean', face_colormap='PiYG')
    points_layer.edge_color_mode = 'colormap'
    napari.run()
