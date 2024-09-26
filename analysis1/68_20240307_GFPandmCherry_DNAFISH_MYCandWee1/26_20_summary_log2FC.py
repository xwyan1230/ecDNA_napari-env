import numpy as np
import pandas as pd

# INPUT PARAMETERS
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"

data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder
hue_order = ['GFP', 'mCherry']
gene = pd.read_csv('%s/gene.txt' % data_dir, na_values=['.'], sep='\t')
ctrl_lst = ['F4']

summary = pd.DataFrame(columns=['sample', 'ctrl', 'limit', 'area_nuclear', 'total_area_ecDNA',
                                'total_area_ratio_ecDNA', 'mean_int_DNAFISH', 'n_ecDNA',
                                'peak_relativer', 'fwhm_relativer', 'center_relativer', 'peak_value_relativer',
                                'peak_absoluter', 'fwhm_absoluter', 'center_absoluter', 'peak_value_absoluter'])


def get_log2FC_phenotype(df, feature, ctrl):
    if ctrl == 'GFP':
        return np.log2((df[feature][0]+0.0001)/(df[feature][1]+0.0001))
    elif ctrl == 'mCherry':
        return np.log2((df[feature][1]+0.0001)/(df[feature][0]+0.0001))


def get_phenotype(df, feature, ctrl):
    if ctrl == 'GFP':
        return df[feature][0]-df[feature][1]
    elif ctrl == 'mCherry':
        return df[feature][1]-df[feature][0]


for hue in hue_order:
    print(hue)
    sample_lst = gene[gene['ctrl'].isin([hue, 'ctrl'])]['sample'].tolist()
    print(sample_lst)
    for sample in sample_lst:
        df = pd.read_csv('%s/%s_summary.txt' % (data_dir, sample), na_values=['.'], sep='\t')
        summary.loc[len(summary.index)] = [sample, hue, df['limit'][2],
                                           get_log2FC_phenotype(df, 'mean_area_nuclear', hue),
                                           get_log2FC_phenotype(df, 'mean_total_area_ecDNA', hue),
                                           get_log2FC_phenotype(df, 'mean_total_area_ratio_ecDNA', hue),
                                           get_log2FC_phenotype(df, 'mean_mean_int_DNAFISH', hue),
                                           get_log2FC_phenotype(df, 'mean_n_ecDNA', hue),
                                           get_log2FC_phenotype(df, 'peak_relativer', hue),
                                           get_log2FC_phenotype(df, 'fwhm_relativer', hue),
                                           get_log2FC_phenotype(df, 'center_relativer', hue),
                                           get_log2FC_phenotype(df, 'peak_value_relativer', hue),
                                           get_log2FC_phenotype(df, 'peak_absoluter', hue),
                                           get_log2FC_phenotype(df, 'fwhm_absoluter', hue),
                                           get_log2FC_phenotype(df, 'center_absoluter', hue),
                                           get_log2FC_phenotype(df, 'peak_value_absoluter', hue)]
summary.to_csv('%ssummary_log2FC.txt' % output_dir, index=False, sep='\t')

for hue in hue_order:
    print(hue)
    summary_ctrl = summary[(summary['ctrl'] == hue) & (summary['sample'].isin(ctrl_lst))].copy().reset_index(drop=True)
    ctrl = pd.DataFrame(columns=['limit', 'area_nuclear', 'total_area_ecDNA', 'total_area_ratio_ecDNA', 'mean_int_DNAFISH',
                                 'n_ecDNA', 'peak_relativer', 'fwhm_relativer', 'center_relativer', 'peak_value_relativer',
                                 'peak_absoluter', 'fwhm_absoluter', 'center_absoluter', 'peak_value_absoluter'])
    ctrl.loc[0] = summary_ctrl.mean(axis=0).tolist()

    ctrl.to_csv('%sctrl_log2FC_%s.txt' % (output_dir, hue), index=False, sep='\t')

print("DONE!")