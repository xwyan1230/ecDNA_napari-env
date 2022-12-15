import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221213_analysis_DNA-FISH_mixing_exp/20221017_mixing-test_mCherry-series_after-heating/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

df = pd.read_csv(("%sDM_mix_DM-H2B-mCherry_ecDNA.txt" % data_dir), na_values=['.'], sep='\t')

hist_colors = [(0.95, 0.50, 0.50), (0.90, 0.90, 0.90), (0.50, 0.90, 0.90)]
line_colors = [(0.85, 0.35, 0.25), (0.30, 0.30, 0.30), (0.30, 0.70, 0.70)]
# '#FFA500', '#40E0D0', '#DA70D6', '#00CED1'
feature = ['percentage_area_curve_ecDNA', 'percentage_area_ratio_curve_ecDNA', 'percentage_int_curve_ecDNA',
           'cum_area_ind_ecDNA', 'cum_area_ratio_ind_ecDNA', 'cum_int_ind_ecDNA', 'cum_area_ind_ecDNA_filled',
           'cum_area_ratio_ind_ecDNA_filled', 'cum_int_ind_ecDNA_filled']

for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

sample = []
for i in range(len(df)):
    if df['mCherry_mean'].tolist()[i] < 5000:
        sample.append('neg')
    elif df['mCherry_mean'].tolist()[i] > 10000:
        sample.append('pos')
    else:
        sample.append('NA')
df['sample'] = sample

df_sort = df[df['sample'].isin(['neg', 'pos'])].copy().reset_index(drop=True)
df_pos = df[df['sample'].isin(['pos'])].copy().reset_index(drop=True)
df_neg = df[df['sample'].isin(['neg'])].copy().reset_index(drop=True)

# cumulative curve
print("Plotting cumulative curve...")
feature = ['cum_int_ind_ecDNA_filled', 'cum_area_ind_ecDNA_filled',
           'cum_area_ratio_ind_ecDNA_filled', 'percentage_area_curve_ecDNA',
           'percentage_int_curve_ecDNA', 'percentage_area_ratio_curve_ecDNA']

for f in feature:
    x_label = 'number of ecDNA hub'

    number_nuclear1 = len(df_pos)
    number_nuclear2 = len(df_neg)

    mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(dat.list_fill_with_last_num(df_pos[f].tolist()))
    mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(dat.list_fill_with_last_num(df_neg[f].tolist()))

    plt.subplots(figsize=(6, 4))
    for i in range(len(df_pos)):
        plt.plot(df_pos[f][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
    for i in range(len(df_neg)):
        plt.plot(df_neg[f][i], alpha=0.05, color=[line_colors[1][j]+0.05 for j in range(len(line_colors[1]))])
    plt.plot(mean_curve1, color=line_colors[0], label='%s, n=%s' % ('pos', number_nuclear1))
    plt.plot(mean_curve2, color=line_colors[1], label='%s, n=%s' % ('neg', number_nuclear2))
    plt.plot(ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
    plt.plot(ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.plot(ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
    plt.xlabel(x_label)
    plt.ylabel(f)
    plt.legend()
    plt.savefig('%s/%s.pdf' % (output_dir, f))
    plt.close()
