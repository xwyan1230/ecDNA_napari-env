import skimage.io as skio
import napari
import cv2
import shared.image as ima
import imutils
import pandas as pd
import matplotlib.pyplot as plt
import shared.display as dis

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'D7'

### 06
dshape_factor = 0.145  # 766nm for Keyence 10x, 360nm for 512x512 at 63x, 58.6 for 3144x3144 at 63x (0.0765), 180 for 1024x1024 at 63x (0.2349), 110 for 2720x2720x at 60x (nikon, 0.145)
interval = 50
align = pd.read_csv('%s/alignment.txt' % output_dir, na_values=['.'], sep='\t')

# img
img_after_hoechst = skio.imread("%s/%s/03_%s_hoechst_final.tif" % (output_dir, sample, sample), plugin="tifffile")
img = imutils.rotate(img_after_hoechst, angle=180)

# img_search
img_hoechst_DNAFISH = skio.imread("%s/%s/05_%s_DNAFISH_hoechst_stitch.tif" % (output_dir, sample, sample), plugin="tifffile")
ori_shape = img_hoechst_DNAFISH.shape
img_hoechst_DNAFISH1 = cv2.resize(img_hoechst_DNAFISH, dsize=(int(ori_shape[1]*dshape_factor), int(ori_shape[0]*dshape_factor)), interpolation=cv2.INTER_AREA)
img_search = img_hoechst_DNAFISH1

# global search
viewer = napari.Viewer()
viewer.add_image(img, blending='additive', colormap='blue', contrast_limits=[0, 40000])
viewer.add_image(img_search, blending='additive', colormap='green', contrast_limits=[0, 20000])
shapes = viewer.add_shapes(name='Shapes', ndim=2)
napari.run()
poly_data = shapes.data[0]
topleft_target, min_ratio = ima.img_search_global(img, img_search, poly_data, interval)

# local view
print(topleft_target)
print(sample)
"""img_cut = img.copy()
img_cut = img_cut[int(topleft_target[1]):int(topleft_target[1]+img_search.shape[0]), int(topleft_target[0]):int(topleft_target[0]+img_search.shape[1])]

viewer = napari.Viewer()
viewer.add_image(img_cut, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_search, blending='additive', colormap='green', contrast_limits=[0, 65535])
shapes = viewer.add_shapes(name='2nd Global check', ndim=2)
napari.run()"""

# global view
# if len(shapes.data) == 0:
align = pd.read_csv('%s/alignment.txt' % output_dir, na_values=['.'], sep='\t')
if len(align[(align['sample'] == sample) & (align['step'] == 'afterFISH-FISH_global')]) == 0:
    align.loc[len(align.index)] = [sample, 'afterFISH-FISH_global', topleft_target[0], topleft_target[1]]
else:
    location = align[(align['sample'] == sample) & (align['step'] == 'afterFISH-FISH_global')].index
    align.loc[location] = [sample, 'afterFISH-FISH_global', topleft_target[0], topleft_target[1]]
align.to_csv('%s/alignment.txt' % output_dir, index=False, sep='\t')

### 07

interval = 5

# local_search
topleft_target, min_ratio = ima.img_search_local(img, img_search, topleft_target, interval)

# local view
img_cut = img.copy()
img_cut = img_cut[int(topleft_target[1]):int(topleft_target[1]+img_search.shape[0]), int(topleft_target[0]):int(topleft_target[0]+img_search.shape[1])]

"""viewer = napari.Viewer()
viewer.add_image(img_cut, blending='additive', colormap='blue', contrast_limits=[0, 40000])
viewer.add_image(img_search, blending='additive', colormap='green', contrast_limits=[0, 20000])
shapes = viewer.add_shapes(name='2nd Local check', ndim=2)
napari.run()"""

# if len(shapes.data) == 0:
viewer = napari.Viewer()
viewer.add_image(img_cut, blending='additive', colormap='blue', contrast_limits=[0, 40000])
viewer.add_image(img_search, blending='additive', colormap='green', contrast_limits=[0, 20000])
plt.imsave("%s%s/07_%s_afterFISH-FISH_local.tiff" % (output_dir, sample, sample), dis.blending(viewer))
viewer.close()

# global view
img_final, img_search_final = ima.img_align_move(img, img_search, [0, 0], topleft_target)

viewer = napari.Viewer()
viewer.add_image(img_final, blending='additive', colormap='blue', contrast_limits=[0, 40000])
viewer.add_image(img_search_final, blending='additive', colormap='green', contrast_limits=[0, 20000])
plt.imsave("%s%s/07_%s_afterFISH-FISH_global.tiff" % (output_dir, sample, sample), dis.blending(viewer))
viewer.close()

align = pd.read_csv('%s/alignment.txt' % output_dir, na_values=['.'], sep='\t')
if len(align[(align['sample'] == sample) & (align['step'] == 'afterFISH-FISH_local')]) == 0:
    align.loc[len(align.index)] = [sample, 'afterFISH-FISH_local', topleft_target[0], topleft_target[1]]
else:
    location = align[(align['sample'] == sample) & (align['step'] == 'afterFISH-FISH_local')].index
    align.loc[location] = [sample, 'afterFISH-FISH_local', topleft_target[0], topleft_target[1]]
align.to_csv('%s/alignment.txt' % output_dir, index=False, sep='\t')

print("DONE!")


