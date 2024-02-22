import skimage.io as skio
import napari
import imutils
import shared.image as ima
import tifffile as tif
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import shared.display as dis
# from napari_animation import Animation
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230711_analysis_H3K4me3_NPC/data/color_img/TPR/fov1_1543.9267693_904.96859135_tif/"

DNAFISH_val = pd.read_csv('%s/DNAFISH_Values.csv' % master_folder, na_values=['.'], sep='\t')
TPR_val = pd.read_csv('%s/TPR_Values.csv' % master_folder, na_values=['.'], sep='\t')

DNAFISH_val['int'] = [float(DNAFISH_val['X,Y'][i].split(',')[1]) for i in range(len(DNAFISH_val))]
TPR_val['int'] = [float(TPR_val['X,Y'][i].split(',')[1]) for i in range(len(TPR_val))]

x = range(len(DNAFISH_val))
plt.subplots(figsize=(12, 9))
plt.plot(x, DNAFISH_val['int'], color='g', label='MYC DNA FISH')
plt.plot(x, TPR_val['int'], color='r', label='TPR IF')
plt.ylabel('Intensity')
plt.legend()
plt.savefig('%s/line_profile.pdf' % master_folder)
plt.close()
