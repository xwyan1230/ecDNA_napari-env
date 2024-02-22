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

df_growth = pd.DataFrame(columns=['sample', '2hr', '6hr', '24hr'])
line_colors = [(255/255, 222/255, 173/255), (255/255, 228/255, 225/255)] * 5 + [(135/255, 205/255, 50/255), (255/255, 127/255, 80/255)]
sns.set_palette(sns.color_palette(line_colors))

for i in range(len(ctrl)):
    df_growth.loc[len(df_growth.index)] = ['DM_%s' % ctrl[i],
                                           df[(df['plate'] == 'DM_FUCCI_2hr') & (df['sample'] == ctrl[i])]['total'].tolist()[0],
                                           df[(df['plate'] == 'DM_FUCCI_6hr') & (df['sample'] == ctrl[i])]['total'].tolist()[0],
                                           df[(df['plate'] == 'DM_FUCCI_24hr') & (df['sample'] == ctrl[i])]['total'].tolist()[0]]
    df_growth.loc[len(df_growth.index)] = ['HSR_%s' % ctrl[i],
        df[(df['plate'] == 'HSR_FUCCI_2hr') & (df['sample'] == ctrl[i])]['total'].tolist()[0],
        df[(df['plate'] == 'HSR_FUCCI_6hr') & (df['sample'] == ctrl[i])]['total'].tolist()[0],
        df[(df['plate'] == 'HSR_FUCCI_24hr') & (df['sample'] == ctrl[i])]['total'].tolist()[0]]

for w in range(len(wells)):
    if wells[w] not in ctrl:
        temp = pd.DataFrame(columns=['sample', '2hr', '6hr', '24hr'])
        temp.loc[len(temp.index)] = ['DM_%s' % wells[w],
                                               df[(df['plate'] == 'DM_FUCCI_2hr') & (df['sample'] == wells[w])]['total'].tolist()[0],
                                               df[(df['plate'] == 'DM_FUCCI_6hr') & (df['sample'] == wells[w])]['total'].tolist()[0],
                                               df[(df['plate'] == 'DM_FUCCI_24hr') & (df['sample'] == wells[w])]['total'].tolist()[0]]
        temp.loc[len(temp.index)] = ['HSR_%s' % wells[w],
                                               df[(df['plate'] == 'HSR_FUCCI_2hr') & (df['sample'] == wells[w])][
                                                   'total'].tolist()[0],
                                               df[(df['plate'] == 'HSR_FUCCI_6hr') & (df['sample'] == wells[w])][
                                                   'total'].tolist()[0],
                                               df[(df['plate'] == 'HSR_FUCCI_24hr') & (df['sample'] == wells[w])][
                                                   'total'].tolist()[0]]
        df_growth_final = pd.concat([df_growth, temp], axis=0).reset_index(drop=True)
        df_growth_final1 = df_growth_final.drop('sample', axis=1)
        df_growth_final1.index = df_growth_final['sample']
        df_growth_final1 = df_growth_final1.transpose()

        fig, ax = plt.subplots(figsize=(10, 9))
        fig.subplots_adjust(right=0.8)
        sns.lineplot(data=df_growth_final1) # annot=True
        plt.legend(loc=(1.04, 0))
        plt.savefig('%s/%s_growth.pdf' % (master_folder, wells[w]))
        plt.show()

