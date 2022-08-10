import numpy as np
import pandas as pd
import random
import math
import scipy.stats as st
import re


"""
# ---------------------------------------------------------------------------------------------------
# FUNCTIONS for PD.DATAFRAME/LIST/UMANAGER/NUMBER
# ---------------------------------------------------------------------------------------------------
List related:

    list_unwrap
        FUNCTION: unwrap one layer of a list with first element per sublist
        SYNTAX:   list_unwrap(lst: list)
    
    str_to_float
        FUNCTION: transform a string into a list of floats
        SYNTAX:   str_to_float(string: str)
    
    find_pos
        FUNCTION: find the position of the first value in linearly increased list that is larger 
                  than or equal to the given value
        SYNTAX:   find_pos(value: int or float, increase_lst: list)
    
    mean_list
        FUNCTION: calculate mean list from a list of lists and its confidence interval, excluding 0 numbers
        SYNTAX:   mean_list(lst: list)
        
    list_exclude_zero
        FUNCTION: exclude zero from list y and corresponding index in x
        SYNTAX:   list_exclude_zero(x: list, y: list)
    
    str_to_list_of_float
        FUNCTION: transform a string into a list of lists of float
        SYNTAX:   str_to_list_of_float(string: str)
    
    list_smooth
        FUNCTION: smooth a list over neighbouring elements
        SYNTAX: list_smooth(lst: list, smooth_factor: int)
    
    list_circle_smooth
        FUNCTION: smooth a list over neighbouring elements, edge will be connected to perform circle smooth
        SYNTAX:   list_circle_smooth(lst: list, smooth_factor: int)
    
    list_peak_center
        FUNCTION: center the maximum number of a list, shift the other elements accordingly
        SYNTAX:   list_peak_center(lst: list)
    
    list_peak_center_with_control
        FUNCTION: center the maximum number of given list lst, shift lst_ctrl based on maximum from lst
        SYNTAX:   list_peak_center_with_control(lst: list, lst_ctrl: list)
    
    list_sum
        FUNCTION: calculate cumulative list based on original list sum[0:i]
        SYNTAX:   list_sum(lst: list)
    
    list_fill_with_num
        FUNCTION: fill given length to the same maximum length with number filled_num
        SYNTAX:   list_fill_with_num(lst: list, filled_num: float)
    
    list_addup_from_df
        FUNCTION: add up all the elements from multiple lists in df[column]
        SYNTAX:   list_addup_from_df(df: pd.DataFrame, column: str)
    
    list_fill_with_last_num
        FUNCTION: fill given length to the same maximum length with last number
        SYNTAX:   list_fill_with_last_num(lst: list)
        
    list_invert
        FUNCTION: invert a list
        SYNTAX:   list_invert(lst:list)
        
Matrix related:
    matrix_pad_with_zeros
        FUNCTION: pad matrix with trailing zeros
        SYNTAX:   matrix_pad_with_zeros(mat: np.array, r1: int, r2: int)

Dataframe related:
    radial_distribution
        FUNCTION: analyze feature_y distribution based on feature_x binning
        SYNTAX:   radial_distribution(df: pd.DataFrame, feature_x, feature_y, interval: float, feature_max: float)
    
    filter_small_from_lst_in_df
        FUNCTION: filter smaller numbers in lists in df[column]
        SYNTAX:   filter_small_from_lst_in_df(df: pd.DataFrame, column: str, threshold: float)

"""


def list_unwrap(lst: list):
    """
    Unwrap one layer of a list with first element per sublist

    Examples:
    input list: [[a1,a2],[b1,b2],[c1,c2],[d1,d2]]
    output list: [a1,b1,c1,d1]
    :param lst: list, the list to be unwrapped
    :return: out: list
    """
    out = list()
    for i in range(len(lst)):
        out.append(lst[i][0])
    return out


def str_to_float(string: str):
    """
    Transform a string into a list of floats

    Examples:
    input string: (24 characters)
    [5.55, 6.53, 7.35, 8.91]
    output list: (4 elements)
    [5.55, 6.53, 7.35, 8.91]
    :param string: str, string to be converted
    :return: out: list
    """
    if len(string) > 2:
        out = [float(i) for i in string[1:-1].split(', ')]
    else:
        out = []
    return out


def find_pos(value: int or float, increase_lst: list):
    """
    Find the position of the first value in linearly increased list that is larger than or equal to the
    given value

    Note: if value > max(increase_lst), which means such position does not exist in the given list. This
        function will return len(increase_lst), which is the last position of the list + 1

    Usage examples:
    1) used to look for bleach frame
       > bleach_frame = find_pos(bleach_time, time_tseries)

    :param value: int or float
    :param increase_lst: list, has to be a linearly increased list
    :return: out: position of the first value in linearly increased list that is larger than or equal to
                the given value, start from 0

    """
    out = len(increase_lst)-1
    i = 0
    while i < len(increase_lst):
        if value <= increase_lst[i]:
            out = i
            break
        else:
            i += 1

    return out


def mean_list(lst: list):
    """
    Calculate mean list from a list of lists and its confidence interval, excluding 0 and nan numbers

    :param lst: list, input list
                [[...], [...], [...], [...],...[...]]
    :return: out: list
    """
    out = []
    out_lower = []
    out_higher = []
    for i in range(len(lst[0])):
        s = []
        for j in range(len(lst)):
            if (lst[j][i] != 0) & (~math.isnan(lst[j][i])):
                s.append(lst[j][i])
        if len(s) != 0:
            out.append(sum(s)*1.0/len(s))
            ci = st.t.interval(alpha=0.95, df=len(s)-1, loc=np.mean(s), scale=st.sem(s))
            out_lower.append(ci[0])
            out_higher.append(ci[1])
        else:
            out.append(0)
            out_lower.append(0)
            out_higher.append(0)

    return out, out_lower, out_higher


def list_exclude_zero(x: list, y: list):
    """
    Exclude zero from list y and corresponding index in x

    :param x: list
    :param y: list
    :return:
    """
    x_out = []
    y_out = []
    for i in range(len(y)):
        if y[i] != 0:
            x_out.append(x[i])
            y_out.append(y[i])

    return x_out, y_out


def matrix_pad_with_zeros(mat: np.array, r1: int, r2: int):
    """
    Pad matrix with trailing zeros

    :param mat: np.array, input matrix
    :param r1: int, expanding length on first axis
    ;param r2: int, expanding length on second axis
    :return:
    """
    l1 = mat.shape[0]
    l2 = mat.shape[1]
    temp1 = np.zeros((l1, r2))
    temp2 = np.zeros((r1, l2 + r2))
    out = np.concatenate((mat, temp1), axis=1)
    out = np.concatenate((out, temp2), axis=0)

    return out


def radial_distribution(df: pd.DataFrame, feature_x, feature_y, interval: float, feature_max: float):
    """
    Analyze feature_y distribution based on feature_x binning

    :param df: pd.DataFamre, including column of feature_x and feature_y
    :param feature_x: string
    :param feature_y: string
    :param interval: float, bin size
    :param feature_max: float, maximum for analysis, number larger than maximum number will be binned in the last bin
    :return:
    """
    out = []
    mean_feature_y = sum(df[feature_y])/len(df)
    for i in np.arange(0, feature_max, interval):
        if feature_max - i <= interval:
            temp = df[df[feature_x] >= i]
        else:
            temp = df[(df[feature_x] >= i) & (df[feature_x] < i+interval)]
        if len(temp) != 0:
            out.append(temp[feature_y].mean()/mean_feature_y)
        else:
            out.append(0)

    return out


def str_to_list_of_float(string: str, n_element: int):
    """
        Transform a string into a list of lists of float

        Examples:
        input string: (28 characters)
        [[5.55, 6.53], [7.35, 8.91]]
        output list: (a list of 2 elements)
        [[5.55, 6.53], [7.35, 8.91]]
        :param string: str, string to be converted
        ;param n_element: int, number of element in the sub list
        :return: out: list
        """
    out = []
    string_split = re.split(r', |[()]', string[1:-1])
    string_split = list(filter(None, string_split))
    string_split = [float(i) for i in string_split]
    for i in range(int(len(string_split)/n_element)):
        out.append([string_split[j+n_element*i] for j in range(n_element)])
    return out


def list_separation(df: pd.DataFrame, key: str):
    """
    Separate 2-element list into two list (mostly for x,y localization)

    :param df: pd.DataFrame
    :param key: str, column name of given feature
    :return:
    """
    x = []
    y = []
    for i in range(len(df)):
        temp = df[key][i]
        x = x + [temp[j][0] for j in range(len(temp))]
        y = y + [temp[j][1] for j in range(len(temp))]

    return x, y


def list_smooth(lst: list, smooth_factor: int):
    """
    Smooth a list over neighbouring elements

    :param lst: list, input list
    :param smooth_factor: int, number of neighbouring elements used to smooth across
    :return:
    """
    out = []
    for i in range(len(lst)-smooth_factor+1):
        temp = lst[i]
        for j in range(smooth_factor-1):
            temp = temp + lst[i+j+1]
        out.append(temp/smooth_factor)

    return out


def list_circle_smooth(lst: list, smooth_factor: int):
    """
    Smooth a list over neighbouring elements, edge will be connected to perform circle smooth

    :param lst: list, input list
    :param smooth_factor: int, number of neighbouring elements used to smooth across
    :return:
    """
    lst_double = list(lst) + list(lst)
    out = []
    for i in range(len(lst)):
        temp = lst[i]
        for j in range(smooth_factor - 1):
            temp = temp + lst_double[i + j + 1]
        out.append(temp / smooth_factor)

    return out


def list_peak_center(lst: list):
    """
    Center the maximum number of a list, shift the other elements accordingly

    :param lst: list, input list
    :return:
    """
    center = int(len(lst)/2)+1
    peak_index = lst.index(np.max(lst))
    lst_temp = lst[peak_index:] + lst[0:peak_index]
    out = lst_temp[center:] + lst_temp[0: center]

    return out


def list_peak_center_with_control(lst: list, lst_ctrl: list):
    """
    Center the maximum number of given list lst, shift lst_ctrl based on maximum from lst

    :param lst: list, input list
    :param lst_ctrl: list, control list
    :return:
    """
    center = int(len(lst)/2)+1
    peak_index = lst.index(np.max(lst))
    lst_temp = lst[peak_index:] + lst[0:peak_index]
    lst_ctrl_temp = lst_ctrl[peak_index:] + lst_ctrl[0:peak_index]
    out = lst_temp[center:] + lst_temp[0:center]
    out_ctrl = lst_ctrl_temp[center:] + lst_ctrl_temp[0:center]

    return out, out_ctrl


def list_sum(lst: list):
    """
    Calculate cumulative list based on original list sum[0:i]

    :param lst: list, input list
    :return:
    """
    out = []
    for i in range(len(lst)+1):
        out.append(sum(lst[0:i]))

    return out


def list_fill_with_num(lst: list, filled_num: float):
    """
    Fill given length to the same maximum length with number filled_num

    :param lst: list, input list
                [[...], [...], [...], [...],...[...]]
    :param filled_num: float, filling number
    :return:
    """
    out = []
    len_lst = [len(lst[i]) for i in range(len(lst))]
    max_len = max(len_lst)
    for i in range(len(lst)):
        if len(lst[i]) < max_len:
            temp = lst[i] + [filled_num]*(max_len-len(lst[i]))
            out.append(temp)
        else:
            out.append(lst[i])

    return out


def list_addup_from_df(df: pd.DataFrame, column: str):
    """
    Add up all the elements from multiple lists in df[column]

    :param df: pd.DataFrame
    :param column: str, column name
    :return:
    """
    out = []
    for i in range(len(df)):
        out = out + df[column][i]

    return out


def filter_small_from_lst_in_df(df: pd.DataFrame, column: str, threshold: float):
    """
    Filter smaller numbers in lists in df[column]

    :param df: pd.DataFrame
    :param column: str, column name
    :param threshold: float, filter threshold
    :return:
    """
    out = []
    for i in range(len(df)):
        temp = [df[column][i][j] for j in range(len(df[column][i])) if df[column][i][j] > threshold]
        out.append(temp)

    return out


def list_fill_with_last_num(lst: list):
    """
    Fill given length to the same maximum length with last number

    :param lst: list, input list
                [[...], [...], [...], [...],...[...]]
    :return:
    """
    out = []
    len_lst = [len(lst[i]) for i in range(len(lst))]
    if len(len_lst) != 0:
        max_len = max(len_lst)
        for i in range(len(lst)):
            if len(lst[i]) < max_len:
                temp = lst[i] + [lst[i][-1]]*(max_len-len(lst[i]))
                out.append(temp)
            else:
                out.append(lst[i])

    return out


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
