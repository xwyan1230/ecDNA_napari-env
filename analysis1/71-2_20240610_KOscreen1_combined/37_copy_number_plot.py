import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'B7'
hueorder = ['mCherry', 'GFP']

data = pd.read_csv('%s/%s/36_%s_copy_number_group.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
data['sum'] = np.log10(data['area_sum'] * (data['DNAFISH_sum'] - data['bg_sum']))
data['merge'] = np.log10(data['area_merge'] * (data['DNAFISH_merge'] - data['bg_merge']))
data['seg_z'] = np.log10(data['area_seg_z'] * (data['DNAFISH_seg_z'] * data['bg_seg_z']))
data = data.sort_values(by='group', ascending=False).copy().reset_index(drop=True)

features = ['sum', 'merge', 'seg_z']
for f in features:
    fig, ax = plt.subplots(figsize=(4, 6))
    fig.subplots_adjust(left=0.4)
    ax = sns.boxplot(data=data, x='group', y=f, hue_order=hueorder)
    lines = ax.get_lines()
    categories = ax.get_xticks()
    plt.close()

    fig, ax = plt.subplots(figsize=(4, 6))
    fig.subplots_adjust(left=0.4)
    sns.swarmplot(data=data, x='group', y=f, hue_order=hueorder, s=2, color=(255 / 255, 140 / 255, 0 / 255))

    nobs = [len(data[data['group'] == hueorder[i]]) for i in range(len(hueorder))]
    nobs = ["%s" % str(x) for x in nobs]
    for cat in categories:
        # every 4th line at the interval of 6 is median line
        # 0 -> p25 1 -> p75 2 -> lower whisker 3 -> upper whisker 4 -> p50 5 -> upper extreme value
        y = round(lines[4 + cat * 6].get_ydata()[0], 2)
        ax.text(cat, y, '%s (%s)' % (y, nobs[cat]), ha='center', va='center', fontweight='bold', size=7, color='white', bbox=dict(facecolor='#445A64'))

    # p = ttest_ind(data_flt[data_flt['group'] == 'mCherry'][f].tolist(), data_flt[data_flt['group'] == 'GFP'][f].tolist())
    # print(p)
    plt.savefig('%s/%s/37_%s_copy_number_%s.pdf' % (output_dir, sample, sample, f))
    # plt.ylim([7, 9])
    # plt.savefig('%s/%s/20_%s_copy_number_fixscale.pdf' % (output_dir, sample, sample))
    plt.show()

print("DONE!")