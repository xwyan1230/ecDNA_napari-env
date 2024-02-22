import os
import skimage.io as skio
import napari
import matplotlib.pyplot as plt
import numpy as np

# parameters (please change)
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/231019_BRCA008_HER2/"
output_directory = '%s/exported_imgs/' % master_folder
threshold_min0 = 0
threshold_max0 = 65535
threshold_min1 = 0
threshold_max1 = 65535
threshold_min2 = 0
threshold_max2 = 65535

imgs = os.listdir(master_folder)
imgs = filter(lambda i: i[-4:] == '.tif', imgs)
imgs = list(set([i.split('_ch')[0] for i in imgs]))
print(len(imgs))  # number of total images requires attention


# necessary functions
def blending(view):
    """
    For exporting napari view images (compress color images together)

    :param view: napari viewer
    :return:
    """
    blended = np.zeros(view.layers[0].data.shape + (4,))
    for layer in view.layers:
        # normalize data by clims
        normalized_data = (layer.data - layer.contrast_limits[0]) / (
                    layer.contrast_limits[1] - layer.contrast_limits[0])
        colormapped_data = layer.colormap.map(normalized_data.flatten())
        colormapped_data = colormapped_data.reshape(normalized_data.shape + (4,))

        blended = blended + colormapped_data

    blended[..., 3] = 1  # set alpha channel to 1
    blended[blended > 1] = 1
    blended[blended < 0] = 0
    return blended


# main
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

for i in range(len(imgs)):
    print('working on %s/%s...' % (i, len(imgs)))
    ch0 = skio.imread("%s/%s_ch00.tif" % (master_folder, imgs[i]), plugin="tifffile")  # default:blue
    ch1 = skio.imread("%s/%s_ch01.tif" % (master_folder, imgs[i]), plugin="tifffile")  # default:green
    ch2 = skio.imread("%s/%s_ch02.tif" % (master_folder, imgs[i]), plugin="tifffile")  # default:red

    viewer = napari.Viewer()
    viewer.add_image(ch0, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(ch1, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(ch2, blending='additive', colormap='red', contrast_limits=[0, 65535])
    plt.imsave('%s/%s_overlay_fullscale.tiff' % (output_directory, imgs[i]), blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(ch0, blending='additive', colormap='blue', contrast_limits=[threshold_min0, threshold_max0])
    viewer.add_image(ch1, blending='additive', colormap='green', contrast_limits=[threshold_min1, threshold_max1])
    viewer.add_image(ch2, blending='additive', colormap='red', contrast_limits=[threshold_min2, threshold_max2])
    plt.imsave('%s/%s_overlay_scaled.tiff' % (output_directory, imgs[i]), blending(viewer))
    viewer.close()

print("DONE!")