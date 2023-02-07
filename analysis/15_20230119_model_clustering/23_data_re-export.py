import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import matplotlib as mpl
import numpy as np
import math
import os
from matplotlib import cm
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230119_model_clustering/"
data_dir = "%stxt/dataset3/0/" % master_folder
output_dir = "%stxt/dataset3/0/" % master_folder

files = [x for x in os.listdir(data_dir)]
if '.DS_Store' in files:
    files.remove('.DS_Store')

for i in files:
    data = pd.read_csv(("%s%s" % (data_dir, i)), na_values=['.'], sep='\t')
    data.insert(1, 'cen_r', [0] * len(data))
    data.to_csv('%s%s' % (output_dir, i), index=False, sep='\t')
