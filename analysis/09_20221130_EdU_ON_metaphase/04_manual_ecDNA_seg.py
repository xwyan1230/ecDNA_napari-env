from skimage.filters import threshold_local, sobel
from skimage.morphology import binary_dilation, binary_erosion, dilation, disk
from skimage.measure import label, regionprops
import shared.objects as obj
from skimage.segmentation import watershed
import shared.segmentation as seg
import skimage.io as skio
from skimage.morphology import extrema
import shared.image as ima
import seaborn as sns
from scipy import ndimage
import napari
import shared.display as dis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tifffile as tif
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221208_analysis_20221130_EdU_ON_metaphase/"
data_dir = '%s20221130_EdU_ON_metaphase/' % master_folder
output_dir = "%sfigures/" % master_folder
sample = 'point5uM_EdU'
start_fov = 31
end_fov = 43

# DATA EXPORT
export_file = '%s%s.txt' % (output_dir, sample)
if os.path.exists(export_file):
    data = pd.read_csv(export_file, na_values=['.'], sep='\t')
else:
    data = pd.DataFrame(columns=['nuclear', 'FOV',
                                 'n', 'EdU_ind_mean_int', 'EdU_ind_area', 'EdU_centroid',
                                 'MYC_ind_mean_int', 'EdU_chr_mean_int'])


# IMAGING ANALYSIS
def fov_to_str(fov):
    if fov < 10:
        out = '00%s' % fov
    else:
        out = '0%s' % fov
    return out


total_fov = end_fov + 1 - start_fov
nuclear = 0
# 2048x2048
for f in range(total_fov):
    fov = f + start_fov
    print(fov)
    file_name = 'Series%s_RAW' % fov_to_str(fov)
    img_hoechst = skio.imread("%s%s_ch00.tif" % (data_dir, file_name), plugin="tifffile")
    img_FISH = skio.imread("%s%s_ch02.tif" % (data_dir, file_name), plugin="tifffile")
    img_EdU = skio.imread("%s%s_ch01.tif" % (data_dir, file_name), plugin="tifffile")
    img_EdU = np.concatenate([np.zeros(shape=[2048, 6]), img_EdU], axis=1)

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_FISH, blending='additive', colormap='green', contrast_limits=[0, img_FISH.max()])
    viewer.add_image(img_EdU, blending='additive', colormap='magenta', contrast_limits=[0, img_EdU.max()])
    shapes = viewer.add_shapes(name='Shapes', ndim=2)
    napari.run()
    shapes_layer = viewer.layers['Shapes']

    # analysis on each separate cell
    for i in range(len(shapes.data)):
        nuclear = nuclear + 1
        poly_data = shapes.data[i]

        # ***** generate local images *****
        # reshape sub_masks
        top, left = np.floor(np.min(poly_data, axis=0))
        bottom, right = np.ceil(np.max(poly_data, axis=0))
        top, bottom = np.clip((top, bottom), 0, img_FISH.shape[0] - 1).astype(int)
        left, right = np.clip((left, right), 0, img_FISH.shape[1] - 1).astype(int)
        output_shape = (bottom - top + 1, right - left + 1)
        # generate sub_masks and sub_channels
        sub_masks = shapes_layer._data_view.to_masks(mask_shape=output_shape, offset=(top, left))[i]
        cell_FISH_mask = img_FISH[top:bottom + 1, left:right + 1] * sub_masks
        cell_hoechst_mask = img_hoechst[top:bottom + 1, left:right + 1] * sub_masks
        cell_EdU_mask = img_EdU[top:bottom + 1, left:right + 1] * sub_masks
        cell_FISH = img_FISH[top:bottom + 1, left:right + 1]
        cell_hoechst = img_hoechst[top:bottom + 1, left:right + 1]
        cell_EdU = img_EdU[top:bottom + 1, left:right + 1]
        # didn't perform background correction since data collected from confocal sp8

        tif.imwrite("%sindividual_cell_tif/%s_%s_%s_hoechst.tif" % (output_dir, sample, fov, nuclear), cell_hoechst)
        tif.imwrite("%sindividual_cell_tif/%s_%s_%s_FISH.tif" % (output_dir, sample, fov, nuclear), cell_FISH)
        tif.imwrite("%sindividual_cell_tif/%s_%s_%s_EdU.tif" % (output_dir, sample, fov, nuclear), cell_EdU)

        # ***** chromosome segmentation *****
        # generate primary chromosome segmentation using EdU channel
        img_for_chr_seg = cell_EdU.copy()
        img_mask_for_chr_seg = cell_EdU_mask.copy()

        local_chr_seg = threshold_local(img_for_chr_seg, 71)  # original 31 for Aarohi
        chromosome_seg = img_mask_for_chr_seg > local_chr_seg
        chromosome_seg = binary_erosion(chromosome_seg)
        chromosome_seg = obj.remove_large(chromosome_seg, 10000)  # original 5000 for Aarohi
        chromosome_seg = binary_erosion(chromosome_seg)
        chromosome_seg = obj.remove_small(chromosome_seg, 100)
        chromosome_seg = ndimage.binary_fill_holes(chromosome_seg)
        chromosome_seg = binary_dilation(chromosome_seg)
        tif.imwrite("%sindividual_cell_tif/%s_%s_%s_chr_seg.tif" % (output_dir, sample, fov, nuclear), chromosome_seg)

        # chromosome seg used to filter FISH signals
        chromosome_seg_filter = binary_dilation(chromosome_seg, disk(5))

        # ***** ecDNA signal segmentation *****
        viewer = napari.Viewer()
        viewer.add_image(cell_FISH, blending='additive', colormap='green', contrast_limits=[0, img_FISH.max()])
        points = viewer.add_points(name='Points', ndim=2)
        napari.run()
        points_lst = points.data

        ecDNA_seg = np.zeros_like(cell_FISH)
        for j in points_lst:
            ecDNA_seg[int(j[0]), int(j[1])] = 1
        ecDNA_seg[chromosome_seg_filter == 1] = 0
        ecDNA_seg[cell_FISH_mask == 0] = 0
        tif.imwrite("%sindividual_cell_tif/%s_%s_%s_ecDNA_seg_by_FISH.tif" % (output_dir, sample, fov, nuclear), ecDNA_seg)

        ecDNA_seg_label = label(ecDNA_seg)
        ecDNA_seg_label = dilation(ecDNA_seg_label, disk(3))

        # MEASUREMENT
        chr_props = regionprops(label(chromosome_seg, connectivity=1), cell_EdU)
        ecDNA_props_EdU = regionprops(ecDNA_seg_label, cell_EdU)
        ecDNA_props_MYC = regionprops(ecDNA_seg_label, cell_FISH)

        EdU_ind_mean_int = [ecDNA_props_EdU[i].intensity_mean for i in range(len(ecDNA_props_EdU))]
        MYC_ind_mean_int = [ecDNA_props_MYC[i].intensity_mean for i in range(len(ecDNA_props_MYC))]

        data.loc[len(data.index)] = [nuclear, fov,
                                     len(ecDNA_props_EdU),
                                     EdU_ind_mean_int,
                                     [ecDNA_props_EdU[i].area for i in range(len(ecDNA_props_EdU))],
                                     [ecDNA_props_EdU[i].centroid for i in range(len(ecDNA_props_EdU))],
                                     MYC_ind_mean_int,
                                     np.mean([chr_props[i].intensity_mean for i in range(len(chr_props))])]
        data.to_csv('%s%s.txt' % (output_dir, sample), index=False, sep='\t')

        EdU_int = pd.DataFrame({'num': range(len(EdU_ind_mean_int)),
                                'EdU_mean': EdU_ind_mean_int,
                                'EdU_normalized_mean': list(np.array(EdU_ind_mean_int)/np.array(MYC_ind_mean_int))})
        # PLOT
        plt.subplots(figsize=(9, 6))
        sns.histplot(data=EdU_int, x='EdU_mean')
        plt.savefig('%shist/%s_%s_%s_EdU_mean.tiff' % (output_dir, sample, fov, nuclear))
        plt.close()

        plt.subplots(figsize=(9, 6))
        sns.histplot(data=EdU_int, x='EdU_normalized_mean')
        plt.savefig('%shist/%s_%s_%s_EdU_normalized_mean.tiff' % (output_dir, sample, fov, nuclear))
        plt.close()

        plt.subplots(figsize=(9, 6))
        plt.scatter(EdU_ind_mean_int, MYC_ind_mean_int)
        plt.savefig('%shist/%s_%s_%s_scatter.tiff' % (output_dir, sample, fov, nuclear))
        plt.show()

        viewer1 = napari.Viewer()
        viewer1.add_image(cell_FISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
        viewer1.add_image(cell_EdU, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
        plt.imsave('%scolor_img/%s_%s_%s_img.tiff' % (output_dir, sample, fov, nuclear), dis.blending(viewer1))
        viewer1.close()

        viewer2 = napari.Viewer()
        viewer2.add_image(chromosome_seg, blending='additive', colormap='magenta', contrast_limits=[0, 1])
        viewer2.add_image(ecDNA_seg_label, blending='additive', colormap='green', contrast_limits=[0, 1])
        plt.imsave('%scolor_img/%s_%s_%s_seg.tiff' % (output_dir, sample, fov, nuclear), dis.blending(viewer2))
        viewer2.close()

        viewer1 = napari.Viewer()
        viewer1.add_image(cell_hoechst, blending='additive', colormap='blue', contrast_limits=[5000, 40000])
        plt.imsave('%scolor_img/%s_%s_%s_hoechst.tiff' % (output_dir, sample, fov, nuclear), dis.blending(viewer1))
        viewer1.close()

        viewer1 = napari.Viewer()
        viewer1.add_image(cell_FISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
        plt.imsave('%scolor_img/%s_%s_%s_FISH.tiff' % (output_dir, sample, fov, nuclear), dis.blending(viewer1))
        viewer1.close()

        viewer1 = napari.Viewer()
        viewer1.add_image(cell_EdU, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
        plt.imsave('%scolor_img/%s_%s_%s_EdU.tiff' % (output_dir, sample, fov, nuclear), dis.blending(viewer1))
        viewer1.close()

        viewer = napari.Viewer()
        viewer.add_image(cell_FISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
        viewer.add_image(cell_EdU, blending='additive', colormap='magenta', contrast_limits=[0, 65535])

        points = [ecDNA_props_EdU[i].centroid for i in range(len(ecDNA_props_EdU))]
        features = {'EdU_mean': EdU_ind_mean_int,
                    'EdU_normalized_mean': list(np.array(EdU_ind_mean_int)/np.array(MYC_ind_mean_int))}

        text = {'string': '{EdU_normalized_mean:.2f}', 'size': 10, 'color': 'white', 'translation': np.array([-5, 5])}
        points_layer = viewer.add_points(points, features=features, text=text, size=5, edge_width=1,
                                         edge_width_is_relative=False, edge_color='EdU_normalized_mean',
                                         edge_colormap='PiYG',
                                         face_color='EdU_normalized_mean', face_colormap='PiYG')
        points_layer.edge_color_mode = 'colormap'
        napari.run()

print("DONE!")







