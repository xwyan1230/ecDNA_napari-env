import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230216_analysis_BRD4-series/"
data_dir = "%sdata/mitosis/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'DM-Ctrl_mix_mCh-POLR3D'
fov = 3
number = 11
file_name = '%s_mitosis%s_RAW' % (sample, fov)
r = 250


img_hoechst = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
img_DNAFISH = skio.imread("%s%s/%s_ch02.tif" % (data_dir, sample, file_name), plugin="tifffile")
img_mCherry = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
img_laminB = skio.imread("%s%s/%s_ch03.tif" % (data_dir, sample, file_name), plugin="tifffile")

print(img_hoechst.shape)

viewer = napari.Viewer()
viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
viewer.add_image(img_mCherry, blending='additive', colormap='red', contrast_limits=[0, img_mCherry.max()])
viewer.add_image(img_laminB, blending='additive', colormap='magenta', contrast_limits=[0, img_laminB.max()])
points = viewer.add_points(name='Points', ndim=2)
# viewer.dims.ndisplay = 3
napari.run()

center = points.data[0]
x1 = int(center[0]-r) if center[0] >= r else int(0)
x2 = int(center[0]+r) if (center[0] + r) < 3144 else int(3144)
y1 = int(center[1]-r) if center[1] >= r else int(0)
y2 = int(center[1]+r) if (center[1] + r) < 3144 else int(3144)

if not os.path.exists("%s%s/mitotic_img/%s_%s/" % (output_dir, sample, sample, number)):
    os.makedirs("%s%s/mitotic_img/%s_%s/" % (output_dir, sample, sample, number))

if len(img_hoechst.shape) > 2:
    for i in range(img_hoechst.shape[0]):
        viewer = napari.Viewer()
        viewer.add_image(img_hoechst[i, x1:x2, y1:y2], blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
        viewer.add_image(img_DNAFISH[i, x1:x2, y1:y2], blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
        viewer.add_image(img_mCherry[i, x1:x2, y1:y2], blending='additive', colormap='red', contrast_limits=[0, img_mCherry.max()])
        plt.imsave('%s%s/mitotic_img/%s_%s/n%s_fov%s_z%s_withmCh.tiff' % (output_dir, sample, sample, number, number, fov, i), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst[i, x1:x2, y1:y2], blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
        viewer.add_image(img_DNAFISH[i, x1:x2, y1:y2], blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
        plt.imsave('%s%s/mitotic_img/%s_%s/n%s_fov%s_z%s.tiff' % (output_dir, sample, sample, number, number, fov, i), dis.blending(viewer))
        viewer.close()
