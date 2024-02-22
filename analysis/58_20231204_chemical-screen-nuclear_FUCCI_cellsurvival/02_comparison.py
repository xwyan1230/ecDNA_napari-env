import pandas as pd
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231204_analysis_FUCCI-screen/"
samples = ['DM_FUCCI_2hr', 'HSR_FUCCI_2hr', 'DM_FUCCI_6hr', 'HSR_FUCCI_6hr', 'DM_FUCCI_24hr', 'HSR_FUCCI_24hr']
ctrl = ['C3', 'C10', 'D6', 'F3', 'F10']
feature_lst = ['per_survival']

df_seq = pd.read_csv('%sseq.txt' % master_folder, na_values=['.'], sep='\t')

for sample in samples:
    df = pd.read_csv('%s/%s/summary_survival.txt' % (master_folder, sample), na_values=['.'], sep='\t')

    for f in range(len(feature_lst)):
        temp = []
        for i in range(len(df_seq)):
            temp.append(df[df['sample'] == df_seq['well'][i]][feature_lst[f]].tolist()[0])
        df_feature = pd.DataFrame(temp)
        df_feature.index = df_seq['compound'].tolist()
        df_feature.columns = [feature_lst[f]]

        fig, ax = plt.subplots(figsize=(9, 14))
        fig.subplots_adjust(left=0.3)
        ax1 = sns.heatmap(df_feature, cbar=0, linewidths=2, vmax=2, vmin=0, square=True, cmap='coolwarm', annot=False, fmt='.2f') # annot=True
        plt.savefig('%s/%s_survival.pdf' % (master_folder, sample))
        plt.close()

