import skimage.io as skio
from skimage.measure import regionprops, label
from skimage.morphology import binary_dilation, disk
import numpy as np
import pandas as pd
import os
import utilities as uti

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/data/"
data_dir1 = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/notes/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/processed/"
batch = 'point5uM_48hr'
plate = 3

# PARAMETERS
exp = 'Colo320_GrayKinase_%s' % batch
otsu_factor = 1.1
circ_threshold = 0.5
min_size = 300
max_size = 2000
threshold = 3000
n_img = 9

skip = pd.read_csv('%s/skip.txt' % data_dir1, na_values=['.'], sep='\t')  # skip broken images
pd_plates = pd.read_csv('%s/plates.txt' % data_dir1, na_values=['.'], sep='\t')
total_wells = pd_plates[(pd_plates['batch'] == batch) & (pd_plates['plate'] == plate)]['total_wells'].tolist()[0]
ks = np.arange(0, total_wells, 1)

samples_lst = ['XY0%s' % (x + 1) for x in range(9)] + ['XY%s' % (x + 10) for x in range(201)]
wells_lst = ['D%s' % (x + 3) for x in range(20)] + ['E%s' % (x + 3) for x in range(20)][::-1] + \
            ['F%s' % (x + 3) for x in range(20)] + ['G%s' % (x + 3) for x in range(20)][::-1] + \
            ['H%s' % (x + 3) for x in range(20)] + ['I%s' % (x + 3) for x in range(20)][::-1] + \
            ['J%s' % (x + 3) for x in range(20)] + ['K%s' % (x + 3) for x in range(20)][::-1] + \
            ['L%s' % (x + 3) for x in range(20)] + ['M%s' % (x + 3) for x in range(20)][::-1] + \
            ['N%s' % (x + 3) for x in range(10)]
samples = samples_lst[:total_wells]
wells = wells_lst[:total_wells]

if batch == '24hr_density':
    densities = [1] * 20 + [2] * 20 + [3] * 20 + [4] * 20 + [6] * 20 + [8] * 20 + \
            [10] * 20 + [12] * 20 + [14] * 20 + [16] * 20 + [-1] * 10  # -1 is control
elif batch == '48hr_density':
    densities = [1] * 20 + [2] * 20 + [3] * 20 + [4] * 20 + [5] * 20 + [6] * 20 + \
                [7] * 20 + [8] * 20 + [9] * 20 + [10] * 20 + [-1] * 10  # -1 is control
else:
    densities = []

for k in ks:
    sample = samples[k]
    well = wells[k]
    print('sample: %s (%s/%s)' % (sample, k+1, total_wells))

    data = pd.DataFrame()
    fov_lst = []
    fov_hoechst_lst = []
    hoechst_lst = []
    h23_lst = []
    azaleaB5_lst = []
    emiRFP670_lst = []

    for i in range(n_img):
        if '%s-%s-%s-%s' % (batch, plate, sample, i+1) not in skip['samples'].tolist():
            file_name = 'Image_%s_0000%s' % (sample, i+1)
            img_hoechst = uti.img_crop(skio.imread("%s%s/%s_%s/%s/%s_CH1.tif" % (data_dir, batch, batch, plate, sample, file_name), plugin="tifffile")[:, :, 2], i)
            img_green = uti.img_crop(skio.imread("%s%s/%s_%s/%s/%s_CH2.tif" % (data_dir, batch, batch, plate, sample, file_name), plugin="tifffile")[:, :, 1], i)
            img_red = uti.img_crop(skio.imread("%s%s/%s_%s/%s/%s_CH3.tif" % (data_dir, batch, batch, plate, sample, file_name), plugin="tifffile")[:, :, 0], i)
            img_farred = uti.img_crop(skio.imread("%s%s/%s_%s/%s/%s_CH4.tif" % (data_dir, batch, batch, plate, sample, file_name), plugin="tifffile")[:, :, 0], i)
            img_nuclear_seg = uti.cell_seg_fluorescent(img_hoechst, otsu_factor=otsu_factor, max_size=max_size,
                                                       maxima_threshold=1.0001,
                                                       min_size=min_size, circ_thresh=circ_threshold,
                                                       threshold=threshold).astype(int)
            nuclear_props = regionprops(img_nuclear_seg)
            centroids = [nuclear_props[j].centroid for j in range(len(nuclear_props))]
            img_seg = np.zeros_like(img_nuclear_seg)
            for j in range(len(centroids)):
                img_seg[int(centroids[j][0])][int(centroids[j][1])] = 1
            img_seg = binary_dilation(img_seg, disk(7))
            img_seg = uti.label_remove_small_large_resort(label(img_seg), 149, 149)
            label_img = label(img_seg)
            blue_props = regionprops(label_img, img_hoechst)
            green_props = regionprops(label_img, img_green)
            red_props = regionprops(label_img, img_red)
            farred_props = regionprops(label_img, img_farred)
            fov_lst = fov_lst + [i] * len(blue_props)
            fov_hoechst_lst = fov_hoechst_lst + [np.sum(img_hoechst) / uti.img_factor(i)] * len(blue_props)
            hoechst_lst = hoechst_lst + [blue_props[j].intensity_mean for j in range(len(blue_props))]
            h23_lst = h23_lst + [green_props[j].intensity_mean for j in range(len(green_props))]
            azaleaB5_lst = azaleaB5_lst + [red_props[j].intensity_mean for j in range(len(red_props))]
            emiRFP670_lst = emiRFP670_lst + [farred_props[j].intensity_mean for j in range(len(farred_props))]

    data['screen'] = [exp] * len(fov_lst)
    data['plate'] = [plate] * len(fov_lst)
    data['sample'] = [sample] * len(fov_lst)
    data['well'] = [well] * len(fov_lst)
    data['fov'] = fov_lst
    data['hoechst'] = hoechst_lst
    data['H2-3'] = h23_lst
    data['AzaleaB5'] = azaleaB5_lst
    data['emiRFP670'] = emiRFP670_lst
    data['fov_hoechst'] = fov_hoechst_lst

    if batch in ['24hr_density', '48hr_density']:
        density = densities[k]
        data['density'] = [density] * len(fov_lst)

    if not os.path.exists("%s%s/%s_%s/txt/" % (output_dir, batch, batch, plate)):
        os.makedirs("%s%s/%s_%s/txt/" % (output_dir, batch, batch, plate))
    data.to_csv('%s/%s/%s_%s/txt/%s.txt' % (output_dir, batch, batch, plate, sample), index=False, sep='\t')

print("DONE!")