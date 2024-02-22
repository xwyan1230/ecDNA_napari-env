import skimage.io as skio
import napari
import imutils
import shared.image as ima
import matplotlib.pyplot as plt
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230724_analysis_GFP-mCherry-test/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'C2'
interval = 50

# img
img_before_GFP = skio.imread("%s/beforeFISH/%s/%s_GFP.tif" % (data_dir, sample, sample), plugin="tifffile")[:, :, 1]
img_before_mCherry = skio.imread("%s/beforeFISH/%s/%s_mCherry.tif" % (data_dir, sample, sample), plugin="tifffile")[:, :, 0]
img_before = img_before_GFP.copy()
img_before[img_before_mCherry > img_before_GFP] = img_before_mCherry[img_before_mCherry > img_before_GFP]
# img_before1 = imutils.rotate(img_before, angle=90)
img = img_before

viewer = napari.Viewer()
viewer.add_image(img, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_before_GFP, blending='additive', colormap='green', contrast_limits=[0, 65535])
viewer.add_image(img_before_mCherry, blending='additive', colormap='red', contrast_limits=[0, 65535])
shapes = viewer.add_shapes(name='Shapes', ndim=2)
napari.run()
shapes_layer = viewer.layers['Shapes']

nuclear = 0
average_GFP_lst = []
average_mCherry_lst = []
for i in range(len(shapes.data)):
    nuclear = nuclear + 1
    poly_data = shapes.data[i]

    # ***** generate local images *****

    # reshape sub_masks
    top, left = np.floor(np.min(poly_data, axis=0))
    bottom, right = np.ceil(np.max(poly_data, axis=0))
    top, bottom = np.clip((top, bottom), 0, img.shape[0] - 1).astype(int)
    left, right = np.clip((left, right), 0, img.shape[1] - 1).astype(int)
    output_shape = (bottom - top + 1, right - left + 1)
    # generate sub_masks and sub_channels
    sub_masks = shapes_layer._data_view.to_masks(mask_shape=output_shape, offset=(top, left))[i]
    img_before_GFP_mask = img_before_GFP[top:bottom + 1, left:right + 1] * sub_masks
    img_before_mCherry_mask = img_before_mCherry[top:bottom + 1, left:right + 1] * sub_masks
    average_GFP = img_before_GFP_mask.sum()/sub_masks.sum()
    average_mCherry = img_before_mCherry_mask.sum() / sub_masks.sum()
    average_GFP_lst.append(np.log(average_GFP))
    average_mCherry_lst.append(np.log(average_mCherry))
    print(average_GFP)
    print(average_mCherry)

print(average_GFP_lst)
print(average_mCherry_lst)

plt.subplots(figsize=(6, 6))
plt.scatter(y=average_GFP_lst, x=average_mCherry_lst, color=(0.30, 0.30, 0.30), s=5, alpha=0.5)
plt.xlim([min(average_mCherry_lst)-0.1, max(average_mCherry_lst)+0.1])
plt.ylim([min(average_GFP_lst)-0.1, max(average_GFP_lst)+0.1])
# plt.savefig('%s/%s/%s_red-green_scatter.pdf' % (output_dir, sample, sample))
plt.show()




