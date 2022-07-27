import skimage.io as skio
import tifffile as tif

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220708_Natasha_ColoDM_interphase/"
sample = 'DMSO1'
raw_folder = '00_raw'
save_folder = '01_raw'
total_fov = 1
total_z_lst = [23]
ch = '01'

# LOAD IMAGE
if total_fov == 1:
    im_stack = skio.imread("%s%s/%s/%s_RAW_ch%s.tif" % (master_folder, sample, raw_folder, sample, ch),
                           plugin="tifffile")
    total_z = total_z_lst[0]
    for i in range(int(im_stack.shape[0]/total_z)):
        im_stack_fov = im_stack[(i * total_z):(i * total_z + total_z)]
        tif.imwrite('%s%s/%s/%s_RAW_ch%s_fov%s.tif' % (master_folder, sample, save_folder, sample, ch, i), im_stack_fov)
else:
    fov_i = 0
    for fov in range(total_fov):
        im_stack = skio.imread("%s%s/%s/%s_%s_RAW_ch%s.tif" % (master_folder, sample, raw_folder, sample, fov+1, ch),
                               plugin="tifffile")
        total_z = total_z_lst[fov]
        for i in range(int(im_stack.shape[0] / total_z)):
            im_stack_fov = im_stack[(i * total_z):(i * total_z + total_z)]
            tif.imwrite('%s%s/%s/%s_RAW_ch%s_fov%s.tif' % (master_folder, sample, save_folder, sample, ch, fov_i+i),
                        im_stack_fov)
        fov_i = fov_i + int(im_stack.shape[0] / total_z)

print("DONE!")

