import shutil
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/grid/" % master_folder
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231016_forKyle_koscreen1_complete/'

sample = 'G11'
temp = os.listdir('%s%s/' % (data_dir, sample))
for i in temp:
    if (i[-3:] == 'tif') & (i[:(4+len(sample))] == '%s_GFP' % sample):
        if not os.path.exists("%s%s_G/" % (output_dir, sample)):
            os.makedirs("%s%s_G/" % (output_dir, sample))
        shutil.copyfile('%s%s/%s' % (data_dir, sample, i), '%s%s_G/%s' % (output_dir, sample, i))
    elif (i[-3:] == 'tif') & (i[:(4+len(sample))] == '%s_mCh' % sample):
        if not os.path.exists("%s%s_R/" % (output_dir, sample)):
            os.makedirs("%s%s_R/" % (output_dir, sample))
        shutil.copyfile('%s%s/%s' % (data_dir, sample, i), '%s%s_R/%s' % (output_dir, sample, i))
print("DONE!")
