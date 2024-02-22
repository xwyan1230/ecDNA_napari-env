import pandas as pd
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import seaborn as sns
import napari
import nd2

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231121_analysis_FUCCI_livecell/"

img_1 = nd2.imread('%s/data/ColoDM_B5_ColoHSR_D8_12conc_4pos.nd2' % master_folder)
img_2 = nd2.imread('%s/data/ColoDM_B5_ColoHSR_D8_12conc_4pos001.nd2' % master_folder)
img = np.concatenate((img_1, img_2), axis=0)
print(img.shape)
# (t, pos, ch, y, x)
pos = 0

img_H2B = img[:, pos, 0, :, :]
img_GFP = img[:, pos, 1, :, :]
img_mCherry = img[:, pos, 2, :, :]

viewer = napari.Viewer()
viewer.add_image(img_H2B, blending='additive', colormap='magenta', contrast_limits=[250, 5000])
viewer.add_image(img_GFP, blending='additive', colormap='green', contrast_limits=[20, 2000])
viewer.add_image(img_mCherry, blending='additive', colormap='red', contrast_limits=[20, 1500])
# shapes_add = viewer.add_shapes(name='add', ndim=2)
napari.run()