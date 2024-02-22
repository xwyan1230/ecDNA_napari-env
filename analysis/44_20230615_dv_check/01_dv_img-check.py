from mrc import DVFile, imread
import napari
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230615_analysis_Ivy_images/"
# sub_folder = '20210903_CHP212_AurB_NMYC'
sub_folder = '20211021_TR14_AurBCy5_CDK4Green_MYCNRed_21'
# sample = 'CHP212_AurBFITC_NMYCRed'
sample = '20211021_TR14_Dual'

total_fov = 11
start_fov = 1
original_spacing = np.array([0.5, 0.1, 0.1])  # my guess

for f in range(total_fov):
    fov = f + start_fov
    if fov < 10:
        file_name = '%s_00%s' % (sample, fov)
    else:
        file_name = '%s_0%s' % (sample, fov)

    img_3d_nuclear = imread('%s%s/%s_R3D_D3D.dv' % (master_folder, sub_folder, file_name))[0, :, :, :]
    img_3d_aurB = imread('%s%s/%s_R3D_D3D.dv' % (master_folder, sub_folder, file_name))[3, :, :, :]
    img_3d_NMYC = imread('%s%s/%s_R3D_D3D.dv' % (master_folder, sub_folder, file_name))[2, :, :, :]
    img_3d_CDK4 = imread('%s%s/%s_R3D_D3D.dv' % (master_folder, sub_folder, file_name))[1, :, :, :]

    viewer = napari.Viewer()
    viewer.add_image(img_3d_nuclear, blending='additive', colormap='blue', name='nuclei', scale=original_spacing)
    viewer.add_image(img_3d_aurB, blending='additive', colormap='red', name='aurB', scale=original_spacing)
    viewer.add_image(img_3d_NMYC, blending='additive', colormap='green', name='NMYC', scale=original_spacing)
    viewer.add_image(img_3d_CDK4, blending='additive', colormap='magenta', name='CDK4', scale=original_spacing)
    napari.run()