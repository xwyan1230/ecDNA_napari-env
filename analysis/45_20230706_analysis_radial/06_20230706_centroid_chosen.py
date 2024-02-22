import skimage.io as skio
import napari
import imutils
import shared.image as ima
import matplotlib.pyplot as plt
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
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

# sample = 'Colo320DM_acidFISH_lamin_3d'
sample = 'Colo320HSR_acidFISH_lamin_3d'
# prefix = '20230614_acidFISH_lamin_ColoDM_DM'
prefix = '20230614_acidFISH_lamin_ColoHSR_HSR'
total_fov = 6
start_fov = 1
srange = 10

df = pd.read_csv('%s%s/%s_centroid.txt' % (data_dir2, sample, sample), na_values=['.'], sep='\t')
df = df[df['area'] < 25000].copy().reset_index(drop=True)

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
    df_f = df[df['FOV'] == fov].copy().reset_index(drop=True)
    while len(df_f) > 0:
        temp_centroid = [df_f['centroid_x'][0], df_f['centroid_y'][0]]
        temp_df = df_f[(df_f['centroid_x'] > (temp_centroid[0]-srange)) & (df_f['centroid_x'] < (temp_centroid[0]+srange)) &
                       (df_f['centroid_y'] > (temp_centroid[1]-srange)) & (df_f['centroid_y'] < (temp_centroid[1]+srange))].copy()
        df_f = df_f.drop(temp_df.index.tolist()).reset_index(drop=True)
        temp_df = temp_df.reset_index(drop=True)
        if len(temp_df) >= 3:
            data.loc[len(data.index)] = temp_df.loc[int(len(temp_df)/2)]

data.to_csv('%s%s/%s_centroid_chosen.txt' % (output_dir, sample, sample), index=False, sep='\t')
print("DONE!")

