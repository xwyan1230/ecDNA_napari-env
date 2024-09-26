import pandas as pd
from aspose.cells import Workbook
import nd2
import napari
import numpy as np
import shared.image as ima
import tifffile as tif

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'D7'
samples = ['D7']
total_fovs = [16]

"""
sample = 'B4'
samples = ['B4_1_9pos', 'B4_2_7pos']
total_fovs = [9, 7]

sample = 'E2'
samples = ['E2_1_3pos', 'E2_2_13pos']
total_fovs = [3, 13]

sample = 'E3'
samples = ['E3_1_12pos', 'E3_2_4pos']
total_fovs = [12, 4]

sample = 'E4'
samples = ['E4_1_14pos', 'E4_2_2pos']
total_fovs = [14, 2]

sample = 'E5'
samples = ['E5_1_12pos', 'E5_2_4pos']
total_fovs = [12, 4]

sample = 'E8'
samples = ['E8_1_1pos', 'E8_2_4pos', 'E8_3_8pos', 'E8_4_3pos']
total_fovs = [1, 4, 8, 3]

sample = 'E9'
samples = ['E9_1_8pos', 'E9_2_8pos']
total_fovs = [8, 8]

sample = 'E11'
samples = ['E11_1_5pos', 'E11_2_4pos', 'E11_3_7pos']
total_fovs = [5, 4, 7]

sample = 'F7'
samples = ['F7_1_1pos', 'F7_2_15pos']
total_fovs = [1, 15]

sample = 'F9'
samples = ['F9_1_8pos', 'F9_2_8pos']
total_fovs = [8, 8]
"""

### 04
for k in range(len(samples)):
    s = samples[k]
    total_fov = total_fovs[k]
    workbook = Workbook("%s/xml/%s.xml" % (data_dir, s))
    workbook.save("%s/%s/04_%s_xml.txt" % (output_dir, sample, s))
    location = pd.read_csv("%s/%s/04_%s_xml.txt" % (output_dir, sample, s), na_values=['.'], sep='\t')
    pd_loc = pd.DataFrame(columns=['fov', 'x', 'y'])
    for i in range(total_fov):
        x = location['Unnamed: %s' % (15 * i + 12)][3]
        y = location['Unnamed: %s' % (15 * i + 14)][3]
        pd_loc.loc[len(pd_loc.index)] = [i, x, y]
    pd_loc.to_csv('%s/%s/04_%s_location.txt' % (output_dir, sample, s), index=False, sep='\t')
print("pd_loc DONE!")

### 05
pixel_size = 300/2720  # uM

pd_loc = pd.DataFrame()
for k in range(len(samples)):
    s = samples[k]
    pd_loc_temp = pd.read_csv("%s/%s/04_%s_location.txt" % (output_dir, sample, s), na_values=['.'], sep='\t')
    pd_loc = pd.concat([pd_loc, pd_loc_temp], axis=0).reset_index(drop=True)

max_x = int((max(pd_loc['x'])-min(pd_loc['x'])+300)/pixel_size)
max_y = int((max(pd_loc['y'])-min(pd_loc['y'])+300)/pixel_size)
min_x1 = min(pd_loc['x'])
max_y1 = max(pd_loc['y'])

img_stitch = np.zeros((max_y, max_x), np.uint16)
print(img_stitch.shape)

for k in range(len(samples)):
    s = samples[k]
    total_fov = total_fovs[k]
    img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, s))
    print(img_stack.shape)
    pd_loc = pd.read_csv("%s/%s/04_%s_location.txt" % (output_dir, sample, s), na_values=['.'], sep='\t')

    for i in range(total_fov):
        img_hoechst = img_stack[i, :, 0, :, :]
        img_hoechst_merge = img_hoechst.max(axis=0)
        delta_x = int((pd_loc['x'][i] - min_x1)/pixel_size)
        delta_y = int((pd_loc['y'][i] - max_y1)/pixel_size)
        print(delta_x)
        print(delta_y)
        img_stitch = ima.image_paste_to(img_stitch, img_hoechst_merge, [-delta_y, delta_x])

"""viewer = napari.Viewer()
viewer.add_image(img_stitch, blending='additive', contrast_limits=[0, 20000])
napari.run()"""

tif.imwrite("%s/%s/05_%s_DNAFISH_hoechst_stitch.tif" % (output_dir, sample, sample), img_stitch)
print("DONE!")





