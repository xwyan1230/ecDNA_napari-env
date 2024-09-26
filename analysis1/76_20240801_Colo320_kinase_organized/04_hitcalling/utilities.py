import numpy as np
import pandas as pd


def get_mean_and_std(data, feature, filter_save1, filter_save0):
    df = data.copy()
    for i in range(len(filter_save1)):
        df = df[df[filter_save1[i]] == 1].copy().reset_index(drop=True)
    for i in range(len(filter_save0)):
        df = df[df[filter_save0[i]] == 0].copy().reset_index(drop=True)
    ctrl_mean = np.mean(df[df['treatment'] == 'DMSO'][feature])
    ctrl_std = np.std(df[df['treatment'] == 'DMSO'][feature])
    return ctrl_mean, ctrl_std


def call_hit(data, feature, filter_save1, filter_save0):
    df = data.copy()
    for i in range(len(filter_save1)):
        df = df[df[filter_save1[i]] == 1].copy()
    for i in range(len(filter_save0)):
        df = df[df[filter_save0[i]] == 0].copy()
    keep_list = df.index
    ctrl_mean, ctrl_std = get_mean_and_std(data, feature, filter_save1, filter_save0)
    print(ctrl_mean)
    print(ctrl_std)
    hit_list = []
    for i in range(len(data)):
        if i not in keep_list:
            hit_list.append(-100)
        elif data[feature][i] < (ctrl_mean - 3*ctrl_std):
            hit_list.append(-1)
        elif data[feature][i] > (ctrl_mean + 3*ctrl_std):
            hit_list.append(1)
        else:
            hit_list.append(0)
    return hit_list


def gene_table(df):
    data = df.copy().reset_index(drop=True)
    gene_lst = []
    for i in range(len(data)):
        gene_lst = gene_lst + list(str(data['target'][i]).split(', '))
    gene_lst = list(set(gene_lst))
    n_gene = []
    for i in range(len(gene_lst)):
        n_gene.append(len(data[data['target'].str.contains(gene_lst[i])]))
    df = pd.DataFrame()
    df['gene'] = gene_lst
    df['n'] = n_gene
    return df


def get_gene_stats(lst, target_df, n_target_limit=5):
    gene_df = gene_table(target_df[target_df['n_target'] < n_target_limit])

    df = pd.DataFrame()
    df['treatment'] = lst
    df['target'] = ['DMSO' if df['treatment'][i] == 'DMSO' else
                    target_df[target_df['treatment'] == df['treatment'][i]]['target'].tolist()[0] for i in
                    range(len(df))]
    df['n_target'] = [len(df['target'][i].split(', ')) for i in range(len(df))]
    df['multi_gene'] = [target_df[target_df['treatment'] == df['treatment'][i]]['multi_gene'].tolist()[0] for i in
                        range(len(df))]
    df_flt = df[df['n_target'] < n_target_limit].copy().reset_index(drop=True)
    df_flt = df_flt[df_flt['multi_gene'] == 1].copy().reset_index(drop=True)
    data = gene_table(df_flt)
    data['n_total'] = [gene_df[gene_df['gene'] == data['gene'][i]]['n'].tolist()[0] for i in range(len(data))]
    data['percentage'] = data['n'] / data['n_total']
    data = data[data['n_total'] > 1].copy().reset_index(drop=True)
    data = data.sort_values(by=['n', 'percentage'], ascending=[False, False]).reset_index(drop=True)
    data_n = data.copy()
    data_n.index = data_n['gene']
    data_n = data_n.drop(['gene', 'n_total', 'percentage'], axis=1)
    data_percentage = data.copy()
    data_percentage.index = data_percentage['gene']
    data_percentage = data_percentage.drop(['gene', 'n_total', 'n'], axis=1)
    return data, data_n, data_percentage


def get_df_t(df, features, batches):
    df_t = pd.DataFrame()
    for batch in batches:
        df_t_temp = pd.DataFrame()
        df_t_temp['group'] = df['group']
        df_t_temp['treatment'] = df['treatment']
        for feature in features:
            df_t_temp[feature] = df['%s_%s' % (batch, feature)]
        df_t_temp['batch'] = [batch] * len(df_t_temp)
        df_t = pd.concat([df_t, df_t_temp], axis=0).reset_index(drop=True)
    return df_t


def get_df_flt(df, filters, filters1, batches):
    df_flt = df.copy()
    if len(filters) != 0:
        for f in filters:
            for batch in batches:
                df_flt = df_flt[df_flt['%s_%s' % (batch, f)] == 1].copy().reset_index(drop=True)
    if len(filters1) != 0:
        for f in filters1:
            for batch in batches:
                df_flt = df_flt[df_flt['%s_%s' % (batch, f)] == 0].copy().reset_index(drop=True)
    return df_flt


def sort_df(df, group_seq):
    df_sort = pd.DataFrame(columns=df.columns)
    for i in range(len(group_seq)):
        df_temp = df[df['group'] == group_seq[i]].copy().reset_index(drop=True)
        if len(df_temp) == 1:
            df_sort.loc[len(df_sort.index)] = df_temp.loc[0]
    return df_sort

