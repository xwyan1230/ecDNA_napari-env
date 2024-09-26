import nd2
import napari
import pandas as pd
import numpy as np
import shared.image as ima
import tifffile as tif
import matplotlib.pyplot as plt
import skimage.io as skio
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240715_analysis_PCC_cellcycle/"
data_dir = '%sprocessed/' % master_folder
output_dir = '%sfigures/' % master_folder

cc = pd.read_excel('%s/cellcycle.xlsx' % master_folder, na_values=['.'])

print(np.mean(cc[cc['group'] == 'M']['n']))
print(np.mean(cc[cc['group'] == 'G1']['n']))
print(np.mean(cc[cc['group'] == 'S']['n']))

print(len(cc[(cc['group'] == 'M')&(cc['n']!=0)])/len(cc[cc['group'] == 'M']))
print(len(cc[(cc['group'] == 'G1')&(cc['n']!=0)])/len(cc[cc['group'] == 'G1']))
print(len(cc[(cc['group'] == 'S')&(cc['n']!=0)])/len(cc[cc['group'] == 'S']))

print(len(cc[cc['group'] == 'M']))
print(len(cc[cc['group'] == 'G1']))
print(len(cc[cc['group'] == 'S']))