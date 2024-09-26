import numpy as np
import pandas as pd


def get_mean_val(data, feature):
    return np.mean(data[(data['treatment'] == 'DMSO') & (data['rep'] == 1)][feature]), np.mean(data[(data['treatment'] == 'DMSO') & (data['rep'] == 2)][feature])


def get_average_delta_log2fc(data_group, data, feature):
    feature_rep1 = data_group[feature][0]
    feature_rep2 = data_group[feature][1]
    feature_mean = (feature_rep1 + feature_rep2)/2
    mean_val_feature = get_mean_val(data, feature)
    delta_feature_rep1 = feature_rep1 - mean_val_feature[0]
    delta_feature_rep2 = feature_rep2 - mean_val_feature[1]
    delta_feature_mean = (delta_feature_rep1 + delta_feature_rep2)/2
    log2fc_feature_rep1 = np.log2((feature_rep1+0.000000001)/(mean_val_feature[0]+0.000000001))
    log2fc_feature_rep2 = np.log2((feature_rep2+0.000000001)/(mean_val_feature[1]+0.000000001))
    log2fc_feature_mean = (log2fc_feature_rep1 + log2fc_feature_rep2)/2
    return feature_rep1, feature_rep2, feature_mean, delta_feature_rep1, delta_feature_rep2, delta_feature_mean, log2fc_feature_rep1, log2fc_feature_rep2, log2fc_feature_mean


def get_cc_percentage(data):
    data['cc_filter'] = [1 if data['check_green'][i] < 50 else 0 for i in range(len(data))]
    data['per_G1'] = data['n_G1']/data['n_filtered']
    data['per_G1S'] = data['n_G1S'] / data['n_filtered']
    data['per_S'] = data['n_S'] / data['n_filtered']
    data['per_G2M'] = data['n_G2M'] / data['n_filtered']
    data['per_G2MG1'] = data['n_G2MG1'] / data['n_filtered']
    data['per_neg_G1'] = data['n_neg_G1'] / data['n_neg']
    data['per_neg_G1S'] = data['n_neg_G1S'] / data['n_neg']
    data['per_neg_S'] = data['n_neg_S'] / data['n_neg']
    data['per_neg_G2M'] = data['n_neg_G2M'] / data['n_neg']
    data['per_neg_G2MG1'] = data['n_neg_G2MG1'] / data['n_neg']
    data['per_pos_G1'] = data['n_pos_G1'] / data['n_pos']
    data['per_pos_G1S'] = data['n_pos_G1S'] / data['n_pos']
    data['per_pos_S'] = data['n_pos_S'] / data['n_pos']
    data['per_pos_G2M'] = data['n_pos_G2M'] / data['n_pos']
    data['per_pos_G2MG1'] = data['n_pos_G2MG1'] / data['n_pos']
    return data


def list_invert(lst: list):
    """
    Invert a list

    :param lst: list, input list
    :return:
    """
    out = []
    lst_temp = lst.copy()
    for i in range(len(lst)):
        out.append(lst_temp[-1])
        lst_temp.pop()

    return out


def get_seq(group_lst, nodes):
    pd_label = pd.DataFrame()
    pd_label['group'] = group_lst
    pd_label_sort = pd.DataFrame(columns=['group'])
    for i in range(len(pd_label)):
        pd_label_sort.loc[len(pd_label_sort.index)] = pd_label.iloc[int(list_invert(nodes)[i])]
    return pd_label_sort['group'].tolist()


def sort_df(df, group_seq):
    df_sort = pd.DataFrame(columns=df.columns)
    for i in range(len(group_seq)):
        df_temp = df[df['group'] == group_seq[i]].copy().reset_index(drop=True)
        if len(df_temp) == 1:
            df_sort.loc[len(df_sort.index)] = df_temp.loc[0]
    return df_sort


def distance(x, y, x0, y0, norm):
    """
    Return distance between point
    P[x0,y0] and a curve (x,y)
    """
    d_x = (x - x0)/norm[0]
    d_y = (y - y0)/norm[1]
    dis = np.sqrt(d_x ** 2 + d_y ** 2)
    return dis


def min_distance(x, y, P, norm, precision=5):
    """
    Compute minimum/a distance/s between
    a point P[x0,y0] and a curve (x,y)
    rounded at `precision`.

    ARGS:
        x, y      (array)
        P         (tuple)
        precision (int)

    Returns min indexes and distances array.
    """
    # compute distance
    d = distance(x, y, P[0], P[1], norm)
    d = np.round(d, precision)
    # find the minima
    glob_min_idxs = np.argwhere(d == np.min(d)).ravel()
    return glob_min_idxs, d