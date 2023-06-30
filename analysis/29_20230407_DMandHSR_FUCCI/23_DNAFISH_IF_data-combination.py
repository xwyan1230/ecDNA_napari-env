import pandas as pd
from shared.sinaplot import sinaplot

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'DM_3_49pos'
df_label = pd.read_csv('%s%s_alignment_label.txt' % (output_dir, sample), na_values=['.'], sep='\t')
df = pd.read_csv('%s%s_n4.txt' % (output_dir, sample), na_values=['.'], sep='\t')
df_cc = pd.read_csv('%sDM_324pos_merge_cellcycle_updated1.txt' % output_dir, na_values=['.'], sep='\t')

df['cellcycle_label'] = df_label['index']
df_sort = df[df['cellcycle_label'] != 0].copy().reset_index(drop=True)

int_green = []
int_red = []
int_MYC = []
int_hoechst = []
area_nuclear = []
int_stdev = []
cellcycle = []

for i in range(len(df_sort)):
    print('%s: %s' % (i, df_sort['cellcycle_label'][i]))
    df_temp = df_cc[df_cc['label'] == df_sort['cellcycle_label'][i]]
    int_green.append(df_temp['mean_int_green'].tolist()[0])
    int_red.append(df_temp['mean_int_red'].tolist()[0])
    int_MYC.append(df_temp['mean_int_MYC'].tolist()[0])
    int_hoechst.append(df_temp['mean_int_nuclear'].tolist()[0])
    area_nuclear.append(df_temp['area_nuclear'].tolist()[0])
    int_stdev.append(df_temp['intensity_stdev'].tolist()[0])
    cellcycle.append(df_temp['cellcycle'].tolist()[0])

df_sort['mean_int_green'] = int_green
df_sort['mean_int_red'] = int_red
df_sort['mean_int_MYC'] = int_MYC
df_sort['mean_int_hoechst'] = int_hoechst
df_sort['area_nuclear_IF'] = area_nuclear
df_sort['intensity_stdev'] = int_stdev
df_sort['cellcycle'] = cellcycle

df_sort.to_csv('%s%s_summary.txt' % (output_dir, sample), index=False, sep='\t')

print("DONE!")