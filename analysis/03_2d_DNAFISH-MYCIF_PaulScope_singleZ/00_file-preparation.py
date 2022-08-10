import os
import shutil

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220726_BRDfamily_screen/"
sample_row_lst = ['B', 'C', 'D', 'E', 'F', 'G']
sample_col_lst = ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11']

for sample_row in sample_row_lst:
    for sample_col in sample_col_lst:
        data_path = "%s%s/%s/" % (master_folder, sample_row, sample_col)
        move_path = "%s%s/%s/01_raw/" % (master_folder, sample_row, sample_col)
        if not os.path.exists(move_path):
            os.makedirs(move_path)
        tifs = [x for x in os.listdir(data_path) if x[-4:] == '.tif']
        for tif in tifs:
            shutil.move('%s%s' % (data_path, tif), move_path)

print("DONE!")
