import pandas as pd
import shared.dataframe as dat

# INPUT PARAMETERS
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240205_analysis_sp8_dWee1andiWee1/"
treatment = 'iWee1_62point5nM_24hr'
samples = ['iWee1_62point5nM_24hr_1']

# NO NEED TO CHANGE
data_dir = "%sprocessed/" % master_folder
output_dir = "%sprocessed/" % master_folder
print(treatment)

df = pd.DataFrame()
for sample in samples:
    df1 = pd.read_csv('%s/txt/%s_n0.txt' % (data_dir, sample), na_values=['.'], sep='\t')
    df = pd.concat([df, df1], axis=0).reset_index(drop=True)

feature_lst = ['DNAFISH_seg_label', 'int_r_to_edge', 'int_relative_r']
for f in feature_lst:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

df_r = pd.DataFrame()
# new approach
data = df
df_r_temp = pd.DataFrame()
seg_lst = []
relative_r_lst = []
absolute_r_lst = []
for i in range(len(data)):
    relative_r_lst = relative_r_lst + data['int_relative_r'][i]
    seg_lst = seg_lst + data['DNAFISH_seg_label'][i]
    absolute_r_lst = absolute_r_lst + data['int_r_to_edge'][i]
df_r_temp['relative_r'] = relative_r_lst
df_r_temp['seg'] = seg_lst
df_r_temp['absolute_r'] = absolute_r_lst
df_r_temp['sample'] = [treatment] * len(df_r_temp)
len_sample = len(df_r_temp[df_r_temp['seg'] == 1.0])
print(len_sample)
# df_r_temp['weights'] = [1.0 / len_sample] * len(df_r_temp)
df_r = pd.concat([df_r, df_r_temp], axis=0)
df_r.to_csv('%s/txt/%s_radial.txt' % (output_dir, treatment), index=False, sep='\t')
print("DONE!")
