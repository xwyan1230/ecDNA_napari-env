import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import random
import shared.image as ima
from skimage.morphology import dilation, disk
import napari

# print(random.choices(np.arange(20, 200, 1), k=10))
r = 75
local_size = 150

im_test = np.zeros((2*local_size, 2*local_size))
im_FISH = np.zeros((2*local_size, 2*local_size))
im_seg = ima.logical_ellipse(im_test, local_size, local_size, r, r)
ones_list = []

for i in range(2*local_size):
    for j in range(2*local_size):
        if im_seg[i][j] == 1:
            ones_list.append([i, j])

random_ones = random.sample(ones_list, 100)
print(random_ones)

for i in random_ones:
    im_FISH[i[0]][i[1]] = 1

im_FISH = dilation(im_FISH, disk(2))
im_FISH[im_seg == 0] = 0

viewer = napari.Viewer()
viewer.add_image(im_seg, colormap='blue', blending='additive')
viewer.add_image(im_FISH, colormap='green', blending='additive', contrast_limits=[0, 1])
napari.run()
