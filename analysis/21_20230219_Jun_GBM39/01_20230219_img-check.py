import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import tifffile as tif
import shared.image as ima

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230219_analysis_Jun_EGFR_RPAs33p_Edu/"
data_dir1 = "%sdata/" % master_folder
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

start_fov = 9
end_fov = 10

sample = 'gbm39ec hu'
total_fov = end_fov - start_fov + 1

for f in range(total_fov):
    fov = start_fov + f
    print(fov)
    file_name = '40x %s-%s' % (sample, fov)
    if (sample == 'gbm39ec hu') & (fov == 9):
        img_EdU = skio.imread("%s%s/C4-%s.tif" % (data_dir1, sample, file_name), plugin="tifffile")[0]
    else:
        img_EdU = skio.imread("%s%s/C4-%s.tif" % (data_dir1, sample, file_name), plugin="tifffile")
    img_hoechst = skio.imread("%s%s/C2-%s.tif" % (data_dir1, sample, file_name), plugin="tifffile")
    img_DNAFISH = skio.imread("%s%s/C1-%s.tif" % (data_dir1, sample, file_name), plugin="tifffile")
    img_RPA = skio.imread("%s%s/C3-%s.tif" % (data_dir1, sample, file_name), plugin="tifffile")
    img_seg = skio.imread("%s%s/seg_tif/%s_nuclear_seg.tif" % (data_dir2, sample, file_name), plugin="tifffile")
    img_ecseg = skio.imread("%s%s/seg_tif/%s_DNAFISH_seg.tif" % (data_dir2, sample, file_name), plugin="tifffile")
    img_RPA_seg = skio.imread("%s%s/seg_tif/%s_RPA_seg.tif" % (data_dir2, sample, file_name), plugin="tifffile")

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    viewer.add_image(img_RPA, blending='additive', colormap='magenta', contrast_limits=[0, img_RPA.max()])
    viewer.add_image(img_EdU, blending='additive', colormap='red', contrast_limits=[0, img_EdU.max()])
    viewer.add_image(img_seg, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_ecseg, blending='additive', colormap='green', contrast_limits=[0, 1])
    viewer.add_image(img_RPA_seg, blending='additive', colormap='magenta', contrast_limits=[0, 1])
    shapes_seg_remove_nuclear = viewer.add_shapes(name='nuclear seg to be removed', ndim=2)
    shapes_seg_add_nuclear = viewer.add_shapes(name='nuclear seg to be added', ndim=2)
    shapes_seg_remove_DNAFISH = viewer.add_shapes(name='DNAFISH seg to be removed', ndim=2)
    shapes_seg_add_DNAFISH = viewer.add_shapes(name='DNAFISH seg to be added', ndim=2)
    shapes_seg_remove_RPA = viewer.add_shapes(name='RPA seg to be removed', ndim=2)
    shapes_seg_add_RPA = viewer.add_shapes(name='RPA seg to be added', ndim=2)
    napari.run()

    img_seg = ima.napari_add_or_remove(shapes_seg_remove_nuclear.data, 'remove', img_seg)
    img_seg = ima.napari_add_or_remove(shapes_seg_add_nuclear.data, 'add', img_seg)
    img_ecseg = ima.napari_add_or_remove(shapes_seg_remove_DNAFISH.data, 'remove', img_ecseg)
    img_ecseg = ima.napari_add_or_remove(shapes_seg_add_DNAFISH.data, 'add', img_ecseg)
    img_RPA_seg = ima.napari_add_or_remove(shapes_seg_remove_RPA.data, 'remove', img_RPA_seg)
    img_RPA_seg = ima.napari_add_or_remove(shapes_seg_add_RPA.data, 'add', img_RPA_seg)

    tif.imwrite('%s%s/seg_tif/%s_nuclear_seg.tif' % (output_dir, sample, file_name), img_seg)
    tif.imwrite('%s%s/seg_tif/%s_DNAFISH_seg.tif' % (output_dir, sample, file_name), img_ecseg)
    tif.imwrite('%s%s/seg_tif/%s_RPA_seg.tif' % (output_dir, sample, file_name), img_RPA_seg)

    viewer = napari.Viewer()
    viewer.add_image(img_seg, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_ecseg, blending='additive', colormap='green', contrast_limits=[0, 1])
    viewer.add_image(img_RPA_seg, blending='additive', colormap='magenta', contrast_limits=[0, 1])
    plt.imsave('%s%s/color_img/%s_seg.tiff' % (output_dir, sample, file_name), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    viewer.add_image(img_RPA, blending='additive', colormap='magenta', contrast_limits=[0, img_RPA.max()])
    viewer.add_image(img_EdU, blending='additive', colormap='red', contrast_limits=[0, img_EdU.max()])
    plt.imsave('%s%s/color_img/%s_img.tiff' % (output_dir, sample, file_name), dis.blending(viewer))
    viewer.close()





