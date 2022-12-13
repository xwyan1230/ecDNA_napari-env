import matplotlib.pyplot as plt
import skimage.io as skio
from skimage.morphology import binary_dilation, binary_erosion, dilation, disk
import napari
import pandas as pd
import shared.dataframe as dat
import shared.display as dis
import seaborn as sns
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221208_analysis_20221130_EdU_ON_metaphase/figures/01_manual/"
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221208_analysis_20221130_EdU_ON_metaphase/figures/01_manual/"
output_dir = data_dir
sample = 'point5uM_EdU'

# load and format data
df = pd.read_csv(("%s%s.txt" % (master_folder, sample)), na_values=['.'], sep='\t')
feature = ['EdU_ind_mean_int', 'EdU_ind_area', 'MYC_ind_mean_int']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

# data analysis
for i in range(len(df)):
    fov = df['FOV'][i]
    nuclear = df['nuclear'][i]
    print(nuclear)

    chr_seg = skio.imread("%sindividual_cell_tif/%s_%s_%s_chr_seg.tif" % (data_dir, sample, fov, nuclear), plugin="tifffile")
    ecDNA_seg = skio.imread("%sindividual_cell_tif/%s_%s_%s_ecDNA_seg_by_FISH.tif" % (data_dir, sample, fov, nuclear), plugin="tifffile")
    ecDNA_seg = dilation(ecDNA_seg, disk(3))

    viewer = napari.Viewer()
    viewer.add_image(chr_seg, blending='additive', colormap='magenta', contrast_limits=[0, 1])
    viewer.add_image(ecDNA_seg, blending='additive', colormap='green', contrast_limits=[0, 1])
    plt.imsave('%scolor_img/%s_%s_%s_seg.tiff' % (output_dir, sample, fov, nuclear), dis.blending(viewer))
    viewer.close()

print("DONE!")