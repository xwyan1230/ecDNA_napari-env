import shutil
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sdata/" % master_folder

treatment = '48hr_density'
save_pos = 4
total_pos = 9

batch1 = os.listdir('%s%s/%s' % (data_dir, treatment, treatment))
batch1 = [k for k in batch1 if '.' not in k] # 59
print(len(batch1))

batch2 = os.listdir('%s%s/%s-1' % (data_dir, treatment, treatment))
batch2 = [k for k in batch2 if '.' not in k]
print(batch2)

for i in range(len(batch2)):
    if batch2[i] == 'XY01':
        source_dir = "%s%s/%s-1/XY01/" % (data_dir, treatment, treatment)
        dest_dir = "%s%s/%s/XY%s/" % (data_dir, treatment, treatment, len(batch1)) if len(batch1) >= 10 else "%s%s/%s/XY0%s/" % (data_dir, treatment, treatment, len(batch1))
        for j in range(total_pos-save_pos):
            for top, dirs, files in os.walk(source_dir):
                search = "Image_XY01_0000%s" % (j + save_pos + 1)
                for k in range(len(files)):
                    if search in files[k]:
                        if len(files[k].split('_')[1]) == 4:
                            new_file_name = files[k][:8]+str(len(batch1))+files[k][10:] \
                                if len(batch1) >= 10 else files[k][:9]+str(len(batch1))+files[k][10:]
                        else:
                            new_file_name = files[k][:8] + str(len(batch1)) + files[k][11:] \
                                if len(batch1) >= 10 else files[k][:9] + str(len(batch1)) + files[k][11:]
                        shutil.move(os.path.join(source_dir, files[k]), os.path.join(dest_dir, new_file_name))
    else:
        source_dir = "%s%s/%s-1/%s/" % (data_dir, treatment, treatment, batch2[i])
        dest_folder_number = int(batch2[i][2:])+len(batch1)-1
        print(dest_folder_number)
        dest_dir = "%s%s/%s/XY%s/" % (data_dir, treatment, treatment, dest_folder_number) if dest_folder_number >= 10 else "%s%s/%s/XY0%s/" % (data_dir, treatment, treatment, dest_folder_number)
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)
        for top, dirs, files in os.walk(source_dir):
            for k in range(len(files)):
                if len(files[k].split('_')[1]) == 4:
                    new_file_name = files[k][:8] + str(dest_folder_number) + files[k][10:] \
                        if dest_folder_number >= 10 else files[k][:9] + str(dest_folder_number) + files[k][10:]
                else:
                    new_file_name = files[k][:8] + str(dest_folder_number) + files[k][11:] \
                        if dest_folder_number >= 10 else files[k][:9] + str(dest_folder_number) + files[k][11:]
                shutil.move(os.path.join(source_dir, files[k]), os.path.join(dest_dir, new_file_name))


