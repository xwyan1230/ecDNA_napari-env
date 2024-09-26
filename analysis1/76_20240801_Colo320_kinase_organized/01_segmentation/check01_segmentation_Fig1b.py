# Fig 1b

import skimage.io as skio
import napari
from skimage.measure import regionprops, label
from skimage.morphology import binary_dilation, disk
import numpy as np
import pandas as pd
import utilities as uti

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/data/"
data_dir1 = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/notes/"
# output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/figures/"
batch = '5uM_48hr'
plate = 2
sample = 'XY77'

# PARAMETERS
otsu_factor = 1.1
circ_threshold = 0.5
min_size = 300
max_size = 2000
threshold = 3000
n_img = 9

skip = pd.read_csv('%s/skip.txt' % data_dir1, na_values=['.'], sep='\t')

for i in range(n_img):
    if '%s_%s_%s' % (batch, sample, i+1) not in skip['samples'].tolist():
        print(i)
        file_name = 'Image_%s_0000%s' % (sample, i+1)
        img_hoechst = uti.img_crop(skio.imread("%s%s/%s_%s/%s/%s_CH1.tif" % (data_dir, batch, batch, plate, sample, file_name), plugin="tifffile")[:, :, 2], i)
        img_green = uti.img_crop(skio.imread("%s%s/%s_%s/%s/%s_CH2.tif" % (data_dir, batch, batch, plate, sample, file_name), plugin="tifffile")[:, :, 1], i)
        img_red = uti.img_crop(skio.imread("%s%s/%s_%s/%s/%s_CH3.tif" % (data_dir, batch, batch, plate, sample, file_name), plugin="tifffile")[:, :, 0], i)
        img_farred = uti.img_crop(skio.imread("%s%s/%s_%s/%s/%s_CH4.tif" % (data_dir, batch, batch, plate, sample, file_name), plugin="tifffile")[:, :, 0], i)
        img_nuclear_seg = uti.cell_seg_fluorescent(img_hoechst, otsu_factor=otsu_factor, max_size=max_size,
                                                   maxima_threshold=1.0001,
                                                   min_size=min_size, circ_thresh=circ_threshold,
                                                   threshold=threshold).astype(int)
        nuclear_props = regionprops(img_nuclear_seg)
        centroids = [nuclear_props[j].centroid for j in range(len(nuclear_props))]
        img_seg = np.zeros_like(img_nuclear_seg)
        for j in range(len(centroids)):
            img_seg[int(centroids[j][0])][int(centroids[j][1])] = 1
        img_seg = binary_dilation(img_seg, disk(7))
        img_seg = uti.label_remove_small_large_resort(label(img_seg), 149, 149)

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 40000])
        viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, 65535])
        viewer.add_image(img_red, blending='additive', colormap='red', contrast_limits=[0, 65535])
        viewer.add_image(img_farred, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
        viewer.add_image(img_seg, blending='additive', contrast_limits=[0, 1])
        """if not os.path.exists("%s%s/%s_%s/seg/%s/" % (output_dir, batch, batch, plate, sample)):
            os.makedirs("%s%s/%s_%s/seg/%s/" % (output_dir, batch, batch, plate, sample))
        plt.imsave("%s%s/%s_%s/seg/%s/seg_%s.tiff" % (output_dir, batch, batch, plate, sample, i + 1), dis.blending(viewer))"""
        napari.run()

print("DONE!")