# Fig 2a

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.optimize import curve_fit
import utilities as uti

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
batch = '5uM_48hr'
feature1 = 'log2fc_hoechst_mean'
feature2 = 'log2fc_n_mean'

# PARAMETERS
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
data = pd.read_csv('%s/01_summary/%s_grouping.txt' % (data_dir, batch), na_values=['.'], sep='\t')
data_filter = 'qc_filter'
data_flt = data[data[data_filter] == 1].copy()
print(len(data_flt))
print(max(data_flt.index))


def model_func(x, a, b, c):
    return a*x**3 + b*x**2 + c*x


popt, _ = curve_fit(model_func, data_flt[feature1].tolist(), data_flt[feature2].tolist())
a, b, c = popt
x = np.arange(np.min(data_flt[feature1])-0.01, np.max(data_flt[feature1])+0.01, 0.001)
y = model_func(np.array(x), a, b, c)
norm = [-0.2, model_func(-0.2, a, b, c)]

"""for i in range(len(data_flt)):
    print(data_flt['treatment'][i])
    P = [data_flt[feature1][i], data_flt[feature2][i]]
    min_idxs, dis = uti.min_distance(x, y, P, norm)
    print(x[min_idxs])
    print(y[min_idxs])
    print(dis[min_idxs])

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(x, y, lw=1)
    for idx in min_idxs:
        ax.plot(
            [P[0], x[idx]],
            [P[1], y[idx]],
            '--', lw=1,
            label=f'distance {dis[idx]:.2f}'
        )
    ax.plot(*P, 'or')
    ax.text(
        P[0], P[1],
        f"  P ({P[0]}, {P[1]})",
        ha='left', va='center',
        fontsize=6
    )
    ax.legend()
    plt.show()"""

dis_lst = []
for i in range(len(data_flt)):
    P = [data_flt[feature1].tolist()[i], data_flt[feature2].tolist()[i]]
    min_idxs, dis = uti.min_distance(x, y, P, norm)
    dis_lst.append(dis[min_idxs[0]])
data_flt['normalized_dis'] = dis_lst

dis_mean = np.mean(data_flt['normalized_dis'])
dis_std = np.std(data_flt['normalized_dis'])
data_flt['cytokinesis_hit'] = [1 if (data_flt[feature1].tolist()[i]>-0.05) & (data_flt['normalized_dis'].tolist()[i] > (dis_mean + 3*dis_std)) else 0 for i in range(len(data_flt))]
data['normalized_dis'] = data_flt['normalized_dis']
data['cytokinesis_hit'] = data_flt['cytokinesis_hit']
data['cytokinesis_hit'] = data['cytokinesis_hit'].fillna(0)

data.to_csv('%s/01_summary/%s_grouping.txt' % (output_dir, batch), index=False, sep='\t')

fig, ax = plt.subplots(figsize=(9, 7))
# fig.subplots_adjust(right=0.8)
sns.scatterplot(data=data_flt, x=feature1, y=feature2, s=20, alpha=1, color='#cccccc')
sns.scatterplot(data=data_flt[data_flt['treatment']=='DMSO'], x=feature1, y=feature2, s=20, alpha=1, color=rainboo_colors[6])
plt.plot(x, y, '--', color=rainboo_colors[1])
data_text = data_flt[data_flt['cytokinesis_hit'] == 1].copy().reset_index(drop=True)
sns.scatterplot(data=data_text, x=feature1, y=feature2, alpha=1, s=20, color=rainboo_colors[2])
for i in range(len(data_text)):
    plt.text(x=data_text[feature1][i] + 0.001, y=data_text[feature2][i], s=data_text['treatment'][i],
             size=6, color=rainboo_colors[2])
plt.show()

print("DONE!")