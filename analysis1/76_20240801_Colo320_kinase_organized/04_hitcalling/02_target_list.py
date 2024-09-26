import pandas as pd
import utilities as uti

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/notes/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/notes/"

# PARAMETERS
target_df = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % data_dir, na_values=['.'])
n_target_limit = 5

df = pd.DataFrame()
df['treatment'] = target_df['compound name']
df['target'] = target_df['TargetGene']
df['n_target'] = [len(df['target'][i].split(', ')) for i in range(len(df))]
print(len(df[df['n_target']>=n_target_limit]))

data = uti.gene_table(df[df['n_target']<n_target_limit])
multi_genelst = data[data['n']>=2]['gene'].tolist()
print(len(multi_genelst))

multi_gene = []
for i in range(len(df)):
    targets = df['target'][i].split(', ')
    count = 0
    for j in targets:
        if (count == 0) & (j in multi_genelst):
            multi_gene.append(1)
            count = 1
    if count == 0:
        multi_gene.append(0)
df['multi_gene'] = multi_gene
print(df)
print(len(df[(df['multi_gene'] == 1) & (df['n_target']<5)]))
df.to_csv('%s/targets.txt' % output_dir, index=False, sep='\t')
print("DONE!")