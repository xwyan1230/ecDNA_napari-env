import pandas as pd
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231204_analysis_FUCCI-screen/"
samples = ['DM_FUCCI_2hr', 'HSR_FUCCI_2hr', 'DM_FUCCI_6hr', 'HSR_FUCCI_6hr', 'DM_FUCCI_24hr', 'HSR_FUCCI_24hr']
ctrl = ['C3', 'C10', 'D6', 'F3', 'F10']

df_seq = pd.read_csv('%sseq.txt' % master_folder, na_values=['.'], sep='\t')

df = pd.DataFrame()
for sample in samples:
    temp = pd.read_csv('%s/%s/summary_survival.txt' % (master_folder, sample), na_values=['.'], sep='\t')
    df = pd.concat([df, temp], axis=0).reset_index(drop=True)

wells = list(set(df['sample'].tolist()))

avg_cell = [np.mean(df[(df['plate'] == 'DM_FUCCI_2hr') & (df['sample'].isin(ctrl))]['total'].tolist()),
            np.mean(df[(df['plate'] == 'DM_FUCCI_6hr') & (df['sample'].isin(ctrl))]['total'].tolist()),
            np.mean(df[(df['plate'] == 'DM_FUCCI_24hr') & (df['sample'].isin(ctrl))]['total'].tolist()),
            np.mean(df[(df['plate'] == 'HSR_FUCCI_2hr') & (df['sample'].isin(ctrl))]['total'].tolist()),
            np.mean(df[(df['plate'] == 'HSR_FUCCI_6hr') & (df['sample'].isin(ctrl))]['total'].tolist()),
            np.mean(df[(df['plate'] == 'HSR_FUCCI_24hr') & (df['sample'].isin(ctrl))]['total'].tolist())]

df_growth = pd.DataFrame(columns=['nor_cell', 'sample', 'treatment', 'x'])

for i in range(len(ctrl)):
    df_growth.loc[len(df_growth.index)] = [df[(df['plate'] == 'DM_FUCCI_2hr') & (df['sample'] == ctrl[i])]['total'].tolist()[0]/avg_cell[0], 'DM_FUCCI_2hr', ctrl[i], 1]
    df_growth.loc[len(df_growth.index)] = [df[(df['plate'] == 'DM_FUCCI_6hr') & (df['sample'] == ctrl[i])]['total'].tolist()[0]/avg_cell[1], 'DM_FUCCI_6hr', ctrl[i], 3]
    df_growth.loc[len(df_growth.index)] = [df[(df['plate'] == 'DM_FUCCI_24hr') & (df['sample'] == ctrl[i])]['total'].tolist()[0]/avg_cell[2], 'DM_FUCCI_24hr', ctrl[i], 5]
    df_growth.loc[len(df_growth.index)] = [df[(df['plate'] == 'HSR_FUCCI_2hr') & (df['sample'] == ctrl[i])]['total'].tolist()[0]/avg_cell[3], 'HSR_FUCCI_2hr', ctrl[i], 1.5]
    df_growth.loc[len(df_growth.index)] = [df[(df['plate'] == 'HSR_FUCCI_6hr') & (df['sample'] == ctrl[i])]['total'].tolist()[0]/avg_cell[4], 'HSR_FUCCI_6hr', ctrl[i], 3.5]
    df_growth.loc[len(df_growth.index)] = [df[(df['plate'] == 'HSR_FUCCI_24hr') & (df['sample'] == ctrl[i])]['total'].tolist()[0]/avg_cell[5], 'HSR_FUCCI_24hr', ctrl[i], 5.5]

for w in range(len(wells)):
    if wells[w] not in ctrl:
        temp = pd.DataFrame(columns=['nor_cell', 'sample', 'treatment', 'x'])
        temp.loc[len(temp.index)] = [df[(df['plate'] == 'DM_FUCCI_2hr') & (df['sample'] == wells[w])]['total'].tolist()[0]/avg_cell[0], 'DM_FUCCI_2hr', wells[w], 1]
        temp.loc[len(temp.index)] = [df[(df['plate'] == 'DM_FUCCI_6hr') & (df['sample'] == wells[w])]['total'].tolist()[0]/avg_cell[1], 'DM_FUCCI_6hr', wells[w], 3]
        temp.loc[len(temp.index)] = [df[(df['plate'] == 'DM_FUCCI_24hr') & (df['sample'] == wells[w])]['total'].tolist()[0]/avg_cell[2], 'DM_FUCCI_24hr', wells[w], 5]
        temp.loc[len(temp.index)] = [df[(df['plate'] == 'HSR_FUCCI_2hr') & (df['sample'] == wells[w])]['total'].tolist()[0]/avg_cell[3], 'HSR_FUCCI_2hr', wells[w], 1.5]
        temp.loc[len(temp.index)] = [df[(df['plate'] == 'HSR_FUCCI_6hr') & (df['sample'] == wells[w])]['total'].tolist()[0]/avg_cell[4], 'HSR_FUCCI_6hr', wells[w], 3.5]
        temp.loc[len(temp.index)] = [df[(df['plate'] == 'HSR_FUCCI_24hr') & (df['sample'] == wells[w])]['total'].tolist()[0]/avg_cell[5], 'HSR_FUCCI_24hr', wells[w], 5.5]

        df_growth_final = pd.concat([df_growth, temp], axis=0).reset_index(drop=True)

        print(wells[w])
        fig, ax = plt.subplots(figsize=(10, 7))
        fig.subplots_adjust(right=0.8)
        line_colors = [(255 / 255, 222 / 255, 173 / 255)] * 6
        sns.set_palette(sns.color_palette(line_colors))
        sns.scatterplot(data=df_growth_final[df_growth_final['treatment'].isin(ctrl)], x='x', y='nor_cell',
                        hue='sample', hue_order=samples)
        line_colors = [(135/255, 205/255, 50/255), (255/255, 127/255, 80/255)] * 3
        sns.set_palette(sns.color_palette(line_colors))
        sns.scatterplot(data=df_growth_final[~df_growth_final['treatment'].isin(ctrl)], x='x', y='nor_cell',
                        hue='sample', hue_order=samples)
        plt.legend(loc=(1.04, 0))
        plt.ylim([0, max([2] + list(df_growth_final['nor_cell']+0.1))])
        plt.savefig('%s/%s_growth.pdf' % (master_folder, wells[w]))
        plt.close()

