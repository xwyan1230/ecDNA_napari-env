import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import shared.dataframe as dat
from scipy.stats import ttest_ind

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'B9'
samples = ['B9']
hueorder =['mCherry', 'GFP']
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red

data_rg = pd.read_csv('%s/%s/19_%s_red_green_group.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
data = pd.DataFrame()
for s in samples:
    df = pd.read_csv('%s/%s/17_%s_cluster.txt' % (output_dir, s, s), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0).reset_index(drop=True)

data = data.sort_values(['sample', 'fov', 'label_nuclear'], ascending=[True, True, True]).copy().reset_index(drop=True)
data_rg = data_rg.sort_values(['sample', 'fov', 'label_mean_int'], ascending=[True, True, True]).copy().reset_index(drop=True)

if (data['label_nuclear'].tolist() == data_rg['label_mean_int'].tolist()) & (data['fov'].tolist() == data_rg['fov'].tolist()):
    data['group'] = data_rg['group']
    data_flt = data[data['group'].isin(['GFP', 'mCherry'])].copy().reset_index(drop=True)
    data_flt = data_flt.sort_values(by='group', ascending=False).copy().reset_index(drop=True)
    feature = ['area_ind_ecDNA', 'area_ratio_ind_ecDNA', 'percentage_area_curve_ecDNA', 'cum_area_ind_ecDNA_filled']
    for f in feature:
        data_flt[f] = [dat.str_to_float(data_flt[f][i]) for i in range(len(data_flt))]
    data_flt1 = data_flt[data_flt['total_area_ratio_ecDNA'] > 0.02].copy().reset_index(drop=True)

    features = ['n_ecDNA', 'total_area_ecDNA', 'total_area_ratio_ecDNA']
    f_range = [[-0.5, 25], [-50, 2000], [-0.01, 0.3]]

    for i in range(len(features)):
        feature = features[i]
        fig, ax = plt.subplots(figsize=(4, 6))
        fig.subplots_adjust(left=0.4)
        ax = sns.boxplot(data=data_flt, x='group', y=feature, hue_order=hueorder)
        lines = ax.get_lines()
        categories = ax.get_xticks()
        plt.close()

        fig, ax = plt.subplots(figsize=(4, 6))
        fig.subplots_adjust(left=0.4)
        sns.swarmplot(data=data_flt, x='group', y=feature, hue_order=hueorder, s=1, color=(255 / 255, 140 / 255, 0 / 255))

        nobs = [len(data_flt[data_flt['group'] == hueorder[i]]) for i in range(len(hueorder))]
        nobs = ["%s" % str(x) for x in nobs]
        for cat in categories:
            # every 4th line at the interval of 6 is median line
            # 0 -> p25 1 -> p75 2 -> lower whisker 3 -> upper whisker 4 -> p50 5 -> upper extreme value
            y = round(lines[4 + cat * 6].get_ydata()[0], 2)
            ax.text(cat, y, '%s (%s)' % (y, nobs[cat]), ha='center', va='center', fontweight='bold', size=7, color='white', bbox=dict(facecolor='#445A64'))

        p = ttest_ind(data_flt[data_flt['group'] == 'mCherry'][feature].tolist(), data_flt[data_flt['group'] == 'GFP'][feature].tolist())
        print(p)
        """y, h, col = data_flt['log10_DNAFISH_total_int_merge'].max() + 0.05, 0.05, 'k'
        plt.plot([0, 0, 1, 1], [y, y + h, y + h, y], lw=1, c=col)
        plt.text((0 + 1) * .5, y + h, "ns", ha='center', va='bottom', color=col)"""
        plt.savefig('%s/%s/21_%s_%s.pdf' % (output_dir, sample, sample, feature))
        plt.ylim(f_range[i])
        plt.savefig('%s/%s/21_%s_%s_fixscale.pdf' % (output_dir, sample, sample, feature))
        plt.show()

    features = ['percentage_area_curve_ecDNA', 'cum_area_ind_ecDNA_filled']
    for f in features:
        x_label = 'number of ecDNA hub'
        df1 = data_flt1[data_flt1['group'] == hueorder[0]].copy().reset_index(drop=True)
        df2 = data_flt1[data_flt1['group'] == hueorder[1]].copy().reset_index(drop=True)
        number_nuclear1 = len(df1)
        number_nuclear2 = len(df2)
        mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(dat.list_fill_with_last_num(df1[f].tolist()))
        mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(dat.list_fill_with_last_num(df2[f].tolist()))

        plt.subplots(figsize=(6, 4))
        for i in range(len(df1)):
            plt.plot(df1[f][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
        for i in range(len(df2)):
            plt.plot(df2[f][i], alpha=0.05, color=[line_colors[1][j] + 0.05 for j in range(len(line_colors[1]))])
        plt.plot(mean_curve1, color=line_colors[0], label='%s, n=%s' % (hueorder[0], number_nuclear1))
        plt.plot(mean_curve2, color=line_colors[1], label='%s, n=%s' % (hueorder[1], number_nuclear2))
        plt.plot(ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
        plt.plot(ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
        plt.plot(ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
        plt.plot(ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
        plt.xlabel(x_label)
        plt.ylabel(f)
        plt.legend()
        plt.savefig('%s/%s/21_%s_%s.pdf' % (output_dir, sample, sample, f))
        plt.xlim([-0.5, 25])
        if f == 'cum_area_ind_ecDNA_filled':
            plt.ylim([-50, 2000])
        plt.savefig('%s/%s/21_%s_%s_fixscale.pdf' % (output_dir, sample, sample, f))
        plt.show()

    features = ['per_AUC']
    f_range = [[0.3, 1]]

    for i in range(len(features)):
        feature = features[i]
        fig, ax = plt.subplots(figsize=(4, 6))
        fig.subplots_adjust(left=0.4)
        ax = sns.boxplot(data=data_flt1, x='group', y=feature, hue_order=hueorder)
        lines = ax.get_lines()
        categories = ax.get_xticks()
        plt.close()

        fig, ax = plt.subplots(figsize=(4, 6))
        fig.subplots_adjust(left=0.4)
        sns.swarmplot(data=data_flt1, x='group', y=feature, hue_order=hueorder, s=1,
                      color=(255 / 255, 140 / 255, 0 / 255))

        nobs = [len(data_flt1[data_flt1['group'] == hueorder[i]]) for i in range(len(hueorder))]
        nobs = ["%s" % str(x) for x in nobs]
        for cat in categories:
            # every 4th line at the interval of 6 is median line
            # 0 -> p25 1 -> p75 2 -> lower whisker 3 -> upper whisker 4 -> p50 5 -> upper extreme value
            y = round(lines[4 + cat * 6].get_ydata()[0], 2)
            ax.text(cat, y, '%s (%s)' % (y, nobs[cat]), ha='center', va='center', fontweight='bold', size=7,
                    color='white', bbox=dict(facecolor='#445A64'))

        p = ttest_ind(data_flt1[data_flt1['group'] == 'mCherry'][feature].tolist(),
                      data_flt1[data_flt1['group'] == 'GFP'][feature].tolist())
        print(p)
        """y, h, col = data_flt['log10_DNAFISH_total_int_merge'].max() + 0.05, 0.05, 'k'
        plt.plot([0, 0, 1, 1], [y, y + h, y + h, y], lw=1, c=col)
        plt.text((0 + 1) * .5, y + h, "ns", ha='center', va='bottom', color=col)"""
        plt.savefig('%s/%s/21_%s_%s.pdf' % (output_dir, sample, sample, feature))
        plt.ylim(f_range[i])
        plt.savefig('%s/%s/21_%s_%s_fixscale.pdf' % (output_dir, sample, sample, feature))
        plt.show()
        """plt.subplots(figsize=(6, 4))
        for i in range(len(df1)):
            plt.plot(df1[f][i], alpha=0.05, color=[line_colors[0][j] + 0.05 for j in range(len(line_colors[0]))])
        plt.plot(mean_curve1, color=line_colors[0], label='%s, n=%s' % (hue_order[0], number_nuclear1))
        plt.plot(ci_lower1, color=line_colors[0], linestyle='--', linewidth=0.5)
        plt.plot(ci_higher1, color=line_colors[0], linestyle='--', linewidth=0.5)
        plt.xlabel(x_label)
        plt.ylabel(f)
        plt.legend()
        plt.savefig('%s/%s/%s_%s_%s.pdf' % (output_dir, sample, sample, f, hue_order[0]))
        plt.close()
    
        plt.subplots(figsize=(6, 4))
        for i in range(len(df2)):
            plt.plot(df2[f][i], alpha=0.05, color=[line_colors[1][j] + 0.05 for j in range(len(line_colors[1]))])
        plt.plot(mean_curve2, color=line_colors[1], label='%s, n=%s' % (hue_order[1], number_nuclear2))
        plt.plot(ci_lower2, color=line_colors[1], linestyle='--', linewidth=0.5)
        plt.plot(ci_higher2, color=line_colors[1], linestyle='--', linewidth=0.5)
        plt.xlabel(x_label)
        plt.ylabel(f)
        plt.legend()
        plt.savefig('%s/%s/%s_%s_%s.pdf' % (output_dir, sample, sample, f, hue_order[1]))
        plt.close()"""

    """features = ['area_ind_ecDNA', 'area_ratio_ind_ecDNA']
    data_m = pd.DataFrame(columns=['sample', 'group'] + features)
    for i in range(len(data_flt)):
        for j in range(data_flt['n_ecDNA'][i]):
            temp = []
            for feature in features:
                temp.append(data_flt[feature][i][j])
            data_m.loc[len(data_m.index)] = [sample, data_flt['group'][i]] + temp

    for feature in features:
        fig, ax = plt.subplots(figsize=(4, 6))
        fig.subplots_adjust(left=0.4)
        ax = sns.boxplot(data=data_m, x='group', y=feature, hue_order=hueorder)
        lines = ax.get_lines()
        categories = ax.get_xticks()
        plt.close()

        fig, ax = plt.subplots(figsize=(4, 6))
        fig.subplots_adjust(left=0.4)
        sns.swarmplot(data=data_m, x='group', y=feature, hue_order=hueorder, s=0.5,
                      color=(255 / 255, 140 / 255, 0 / 255))

        nobs = [len(data_m[data_m['group'] == hueorder[i]]) for i in range(len(hueorder))]
        nobs = ["%s" % str(x) for x in nobs]
        for cat in categories:
            y = round(lines[4 + cat * 6].get_ydata()[0], 2)
            ax.text(cat, y, '%s (%s)' % (y, nobs[cat]), ha='center', va='center', fontweight='bold', size=7,
                    color='white', bbox=dict(facecolor='#445A64'))

        plt.savefig('%s/%s/21_%s_%s.pdf' % (output_dir, sample, sample, feature))
        plt.show()"""
    data.to_csv('%s/%s/21_%s_cluster_group.txt' % (output_dir, sample, sample), index=False, sep='\t')

else:
    print('no')

print("DONE!")