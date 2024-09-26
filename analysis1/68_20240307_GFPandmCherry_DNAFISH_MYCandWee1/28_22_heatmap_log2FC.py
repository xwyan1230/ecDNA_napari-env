import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# INPUT PARAMETERS
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"

data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

feature_lst = ['area_nuclear', 'total_area_ecDNA', 'total_area_ratio_ecDNA', 'mean_int_DNAFISH', 'n_ecDNA',
                                'peak_relativer', 'fwhm_relativer', 'center_relativer', 'peak_value_relativer',
                                'peak_absoluter', 'fwhm_absoluter', 'center_absoluter', 'peak_value_absoluter']

df = pd.read_csv('%ssummary_log2FC.txt' % data_dir, na_values=['.'], sep='\t')
df_ctrl_GFP = pd.read_csv('%sctrl_log2FC_GFP.txt' % data_dir, na_values=['.'], sep='\t')
df_ctrl_mCherry = pd.read_csv('%sctrl_log2FC_mCherry.txt' % data_dir, na_values=['.'], sep='\t')
df_seq = pd.read_csv('%sseq.txt' % data_dir, na_values=['.'], sep='\t')
hue_order = ['GFP', 'mCherry']


for hue in hue_order:
    for feature in feature_lst:
        df_hue = df[df['ctrl'] == hue].copy().reset_index(drop=True)
        df_seq_hue = df_seq[df_seq['ctrl'].isin(['ctrl', hue])].copy().reset_index(drop=True)
        temp = []
        for i in range(len(df_seq_hue)):
            temp.append(df_hue[df_hue['sample'] == df_seq_hue['well'][i]][feature].tolist()[0])
        df_feature = pd.DataFrame(temp)
        df_feature.index = df_seq_hue['gene'].tolist()
        df_feature.columns = [feature]

        if hue == 'GFP':
            ctrl = df_ctrl_GFP
        elif hue == 'mCherry':
            ctrl = df_ctrl_mCherry
        feature_mean = ctrl[feature][0]
        feature_delta_max = max(abs(np.max(temp)-feature_mean), abs(feature_mean-np.min(temp)))

        fig, ax = plt.subplots(figsize=(9, 14))
        fig.subplots_adjust(left=0.3)
        ax1 = sns.heatmap(df_feature, cbar=0, linewidths=2, vmax=feature_mean+feature_delta_max+0.1, vmin=feature_mean-feature_delta_max-0.1, square=True, cmap='coolwarm', annot=False, fmt='.2f') # annot=True
        if not os.path.exists("%s/heatmap_line/" % output_dir):
            os.makedirs("%s/heatmap_line/" % output_dir)
        plt.savefig('%s/heatmap_line/%s_log2FC_%s.pdf' % (output_dir, feature, hue))
        plt.close()