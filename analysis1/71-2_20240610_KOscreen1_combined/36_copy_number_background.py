import pandas as pd
import nd2
import skimage.io as skio
import napari
import numpy as np
import shared.image as ima
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'B7'

data_rg = pd.read_csv('%s/%s/19_%s_red_green_group.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
samples = list(set(data_rg['sample']))
n_sample = len(samples)
print(n_sample)
n_mCherry = 0
n_GFP = 0
n_fov = 0

data = pd.DataFrame(columns=['sample', 'fov', 'label_mean_int', 'group', 'area_sum', 'DNAFISH_sum', 'bg_sum',
                             'area_merge', 'DNAFISH_merge', 'bg_merge', 'area_seg_z', 'DNAFISH_seg_z', 'bg_seg_z'])
pd_seg_z = pd.read_csv('%s/%s/10_%s_seg_z.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')

for k in range(len(samples)):
    s = samples[k]
    img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, s))
    pd_seg_z_s = pd_seg_z[pd_seg_z['sample'] == s].copy().reset_index(drop=True)
    data_s = data_rg[data_rg['sample'] == s].copy().reset_index(drop=True)
    fovs = set(data_s['fov'])
    total_fov = len(fovs)
    for fov in fovs:
        print("%s/%s" % (fov + 1, total_fov))
        n_fov = n_fov + 1
        img_DNAFISH = img_stack[fov, :, 1, :, :]
        img_DNAFISH_merge = img_DNAFISH.max(axis=0)
        img_DNAFISH_sum = img_DNAFISH.astype(int).sum(axis=0)
        seg_z = pd_seg_z_s['seg_z'][fov]
        img_DNAFISH_seg_z = img_DNAFISH[seg_z]
        img_seg_total = skio.imread("%s/%s/13_seg_new_tif/%s_%s_seg_new.tif" % (output_dir, sample, s, fov), plugin="tifffile")
        img_seg_z = skio.imread("%s/%s/13_seg_new_tif/%s_%s_seg_z_new.tif" % (output_dir, sample, s, fov), plugin="tifffile")
        data_fov = data_s[data_s['fov'] == fov].copy().reset_index(drop=True)
        data_fov_mCherry = data_fov[data_fov['group'] == 'mCherry'].copy().reset_index(drop=True)
        data_fov_GFP = data_fov[data_fov['group'] == 'GFP'].copy().reset_index(drop=True)

        if len(data_fov_mCherry) != 0:
            data_fov_mCherry = data_fov_mCherry.sample(n=np.min([7*n_fov-n_mCherry, len(data_fov_mCherry)]), random_state=1)
            n_mCherry = n_mCherry + np.min([7*n_fov-n_mCherry, len(data_fov_mCherry)])

        if len(data_fov_GFP) != 0:
            data_fov_GFP = data_fov_GFP.sample(n=np.min([7*n_fov-n_GFP, len(data_fov_GFP)]), random_state=1)
            n_GFP = n_GFP + np.min([7*n_fov-n_GFP, len(data_fov_GFP)])

        data_measure = pd.concat([data_fov_mCherry, data_fov_GFP], axis=0).reset_index(drop=True)

        for i in range(len(data_measure)):
            print("%s/%s" % (i + 1, len(data_measure)))
            label_int = data_measure['label_mean_int'][i]
            img_cell = img_DNAFISH_merge.copy()
            img_cell[img_seg_z != label_int] = 0
            img_cell_seg_total = img_seg_total.copy()
            img_cell_seg_total[img_seg_total != label_int] = 0
            img_cell_seg_z = img_seg_z.copy()
            img_cell_seg_z[img_seg_z != label_int] = 0
            img_seg_bg = np.zeros_like(img_DNAFISH_merge)

            viewer = napari.Viewer()
            viewer.add_image(img_cell, blending='additive', colormap='green', contrast_limits=[0, 10000])
            shapes = viewer.add_shapes(name='Shapes', ndim=2)
            napari.run()

            img_seg_bg = ima.napari_add_or_remove_obj(shapes.data, 'add', img_seg_bg)

            sum_props = regionprops(label(img_cell_seg_total), img_DNAFISH_sum)
            merge_props = regionprops(label(img_cell_seg_total), img_DNAFISH_merge)
            z_props = regionprops(label(img_cell_seg_z), img_DNAFISH_seg_z)

            sum_bg_props = regionprops(label(img_seg_bg), img_DNAFISH_sum)
            merge_bg_props = regionprops(label(img_seg_bg), img_DNAFISH_merge)
            z_bg_props = regionprops(label(img_seg_bg), img_DNAFISH_seg_z)

            data.loc[len(data.index)] = [s, fov, label_int, data_measure['group'][i], sum_props[0].area,
                                         sum_props[0].intensity_mean, sum_bg_props[0].intensity_mean,
                                         merge_props[0].area, merge_props[0].intensity_mean,
                                         merge_bg_props[0].intensity_mean, z_props[0].area, z_props[0].intensity_mean,
                                         z_bg_props[0].intensity_mean]

        print(n_fov)
        print(n_mCherry)
        print(n_GFP)

data.to_csv('%s/%s/36_%s_copy_number_group.txt' % (output_dir, sample, sample), index=False, sep='\t')
print("DONE!")



