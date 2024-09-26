import numpy as np
import matplotlib.pyplot as plt


def cell_cutoff(data, hc):
    data['log10_emiRFP670'] = np.log10(data['emiRFP670'])
    plt.subplots(figsize=(9, 7))
    m = plt.hist(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])]['log10_emiRFP670'],
                 weights=np.ones(len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])])) / len(
                     data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])]),
                 range=[2.5, 5], bins=60, color='w', edgecolor='black')
    test_index = m[0].tolist().index(np.max(m[0][:15]))
    test_val = m[1][test_index]
    cutoff = test_val + 0.3
    plt.close()
    return cutoff


def cc_cutoff(data, hc):
    data['log2_H2-3'] = np.log2(data['H2-3'])
    data['log2_AzaleaB5'] = np.log2(data['AzaleaB5'])
    data = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])].reset_index(drop=True)
    plt.subplots(figsize=(9, 7))
    m = plt.hist(data['log2_H2-3'], bins=50)
    red_index = m[0].tolist().index(np.max(m[0][:-2]))
    red_axis = m[1][red_index]
    if red_index >= 15:
        red_axis = 10.1
    plt.close()

    plt.subplots(figsize=(9, 7))
    m = plt.hist(data['log2_AzaleaB5'], bins=50)
    green_index = m[0].tolist().index(np.max(m[0][:-2]))
    green_axis = m[1][green_index]
    if green_index >= 15:
        green_axis = 8.3
    plt.close()

    green_cutoff = green_axis + 0.65
    red_cutoff = red_axis + 0.65
    check_green = len(data[data['log2_AzaleaB5'] < (green_axis - 0.65)])
    check_red = len(data[data['log2_H2-3'] < (red_axis - 0.65)])
    return red_cutoff, green_cutoff, check_red, check_green


def get_cc(data, red_cutoff, green_cutoff):
    data['log2_H2-3'] = np.log2(data['H2-3'])
    data['log2_AzaleaB5'] = np.log2(data['AzaleaB5'])
    k = (16 - red_cutoff) / (16 - green_cutoff)
    b = (red_cutoff - green_cutoff) * 16 / (16 - green_cutoff)
    data_G1 = data[(data['log2_AzaleaB5'] > green_cutoff) & (data['log2_H2-3'] <= red_cutoff)].copy().reset_index(drop=True)
    data_G1S = data[(data['log2_AzaleaB5'] <= green_cutoff) & (data['log2_H2-3'] < red_cutoff)].copy().reset_index(drop=True)
    data_S = data[(data['log2_AzaleaB5'] < green_cutoff) & (data['log2_H2-3'] >= red_cutoff)].copy().reset_index(drop=True)
    data_G2M = data[(data['log2_AzaleaB5'] >= green_cutoff) & (data['log2_H2-3'] - k * data['log2_AzaleaB5'] - b > 0)].copy().reset_index(drop=True)
    data_G2MG1 = data[(data['log2_H2-3'] > red_cutoff) & (data['log2_H2-3'] - k * data['log2_AzaleaB5'] - b <= 0)].copy().reset_index(drop=True)
    return data_G1, data_G1S, data_S, data_G2M, data_G2MG1


def get_normalization(data, ctrl, ref):
    mean_hoechst = np.mean(ref[ref['n_filtered'] > 1000]['log10_fov_hoechst'])
    mean_n = np.mean(ref[ref['n_filtered'] > 1000]['n_filtered'])
    hoechst_ratio = np.mean(ctrl[ctrl['n_filtered'] > 1000]['log10_fov_hoechst']) / mean_hoechst
    n_ratio = np.mean(ctrl[ctrl['n_filtered'] > 1000]['n_filtered']) / mean_n
    data['log10_fov_hoechst_normalized'] = data['log10_fov_hoechst'] / hoechst_ratio
    data['n_filtered_normalized'] = data['n_filtered'] / n_ratio
    data['n_pos_normalized'] = data['n_pos'] / n_ratio
    data['n_neg_normalized'] = data['n_neg'] / n_ratio
    mean_G1_hoechst = np.mean(ref[ref['n_filtered'] > 1000]['hoechst_G1'])
    G1_ratio = np.mean(ctrl[ctrl['n_filtered'] > 1000]['hoechst_G1']) / mean_G1_hoechst
    mean_G1S_hoechst = np.mean(ref[ref['n_filtered'] > 1000]['hoechst_G1S'])
    G1S_ratio = np.mean(ctrl[ctrl['n_filtered'] > 1000]['hoechst_G1S']) / mean_G1S_hoechst
    mean_S_hoechst = np.mean(ref[ref['n_filtered'] > 1000]['hoechst_S'])
    S_ratio = np.mean(ctrl[ctrl['n_filtered'] > 1000]['hoechst_S']) / mean_S_hoechst
    mean_G2M_hoechst = np.mean(ref[ref['n_filtered'] > 1000]['hoechst_G2M'])
    G2M_ratio = np.mean(ctrl[ctrl['n_filtered'] > 1000]['hoechst_G2M']) / mean_G2M_hoechst
    mean_G2MG1_hoechst = np.mean(ref[ref['n_filtered'] > 1000]['hoechst_G2MG1'])
    G2MG1_ratio = np.mean(ctrl[ctrl['n_filtered'] > 1000]['hoechst_G2MG1']) / mean_G2MG1_hoechst
    data['hoechst_G1_normalized'] = data['hoechst_G1']/G1_ratio
    data['hoechst_G1S_normalized'] = data['hoechst_G1S'] / G1S_ratio
    data['hoechst_S_normalized'] = data['hoechst_S'] / S_ratio
    data['hoechst_G2M_normalized'] = data['hoechst_G2M'] / G2M_ratio
    data['hoechst_G2MG1_normalized'] = data['hoechst_G2MG1'] / G2MG1_ratio
    data['hoechst_neg_G1_normalized'] = data['hoechst_neg_G1'] / G1_ratio
    data['hoechst_neg_G1S_normalized'] = data['hoechst_neg_G1S'] / G1S_ratio
    data['hoechst_neg_S_normalized'] = data['hoechst_neg_S'] / S_ratio
    data['hoechst_neg_G2M_normalized'] = data['hoechst_neg_G2M'] / G2M_ratio
    data['hoechst_neg_G2MG1_normalized'] = data['hoechst_neg_G2MG1'] / G2MG1_ratio
    data['hoechst_pos_G1_normalized'] = data['hoechst_pos_G1'] / G1_ratio
    data['hoechst_pos_G1S_normalized'] = data['hoechst_pos_G1S'] / G1S_ratio
    data['hoechst_pos_S_normalized'] = data['hoechst_pos_S'] / S_ratio
    data['hoechst_pos_G2M_normalized'] = data['hoechst_pos_G2M'] / G2M_ratio
    data['hoechst_pos_G2MG1_normalized'] = data['hoechst_pos_G2MG1'] / G2MG1_ratio
    return data
