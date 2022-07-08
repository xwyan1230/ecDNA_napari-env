import numpy as np
from matplotlib import cm
from vispy.color import Colormap
from matplotlib.colors import ListedColormap
from matplotlib import pyplot as plt
import pandas as pd
from scipy.stats import ks_2samp
import seaborn as sns
import matplotlib.colors as colors

"""
# ---------------------------------------------------------------------------------------------------
# FUNCTIONS for OUTPUT DISPLAY
# ---------------------------------------------------------------------------------------------------

plot_minus_ln_p
    FUNCTION: generate -ln(p) plot for given feature
    SYNTAX:   plot_minus_ln_p(inc: int, limit: int, repeat: int, feature: str, data_pd: pd.DataFrame, ctrl_lst: 
              list, sample: str, save_path: str)

get_p
    FUNCTION: calculate pair wise KS test p-value for -ln(p) plot
    SYNTAX:   get_p(data1: pd.DataFrame, data2: pd.DataFrame, feature: str, inc: int, limit: int, repeat: int)
    
get_x
    FUNCTION: create pair-wise x value for -ln(p) plot
    SYNTAX:   get_x(inc: int, limit: int, repeat: int, offset: float)

plot_volcano_hit
    FUNCTION: plot volcano plot for FRAP screen (with hit displayed)
    SYNTAX:   plot_volcano_hit(pd_p: pd.DataFrame, pd_value: pd.DataFrame, pd_gamma: pd.DataFrame, feature: str, 
              center: float, threshold: float, show_gene: str, save_path: str)
"""


def plot_minus_ln_p(inc: int, limit: int, repeat: int, feature: str, data_pd: pd.DataFrame, ctrl_lst: list,
                    sample: str, save_path: str):
    """
    Generate -ln(p) plot for given feature

    :param inc: int, increment (generally 5)
    :param limit: int, upper limit of the plot
    :param repeat: int, how many runs to be calculated per condition (generally 50)
    :param feature: str, comparing feature, column name
    :param data_pd: pd.DataFrame, total sample
    :param ctrl_lst: list, list of control wells
    :param sample: str, sample name
    :param save_path: str, saving path
    :return:
    """
    sample_lst = ctrl_lst + [sample, 'WT']
    for i in sample_lst:
        n_curve = len(data_pd[data_pd['sample'] == i])
        if n_curve < limit:
            limit = n_curve

    x = np.arange(limit+5)
    plt.figure(figsize=(15, 4))
    n_middle = (len(ctrl_lst) + 1) // 2
    pd_WT = data_pd[data_pd['sample'] == 'WT']
    for i in range(len(ctrl_lst)):
        pd_ctrl = data_pd[data_pd['sample'] == ctrl_lst[i]]
        plt.scatter(get_x(inc, limit, repeat, (i-n_middle)/2.0),
                    get_p(pd_ctrl, pd_WT, feature, inc, limit, repeat),
                    alpha=0.5, s=5, c='#40E0D0', label=ctrl_lst[i])

    pd_sample = data_pd[data_pd['sample'] == sample]
    plt.scatter(get_x(inc, limit, repeat, (len(ctrl_lst)-n_middle)/2.0),
                get_p(pd_sample, pd_WT, feature, inc, limit, repeat),
                alpha=0.5, s=5, c='#FF4500', label=sample)
    plt.plot(x, 0 * x + 5, linestyle='--', color='#696969')
    plt.xlabel('number of samples')
    plt.ylabel('%s -ln(KS p-value)' % feature)
    plt.savefig('%s%s_%s_p.pdf' % (save_path, sample, feature))
    plt.close()


def get_p(data1: pd.DataFrame, data2: pd.DataFrame, feature: str, inc: int, limit: int, repeat: int):
    """
    Calculate pair wise KS test p-value for -ln(p) plot

    :param data1: pd.DataFrame, data1
    :param data2: pd.DataFrame, data2
    :param feature: str, comparing feature, column name
    :param inc: int, increment (generally 5)
    :param limit: int, upper limit of the plot
    :param repeat: int, how many runs to be calculated per condition (generally 50)
    :return: out: list, list of p-value
    """
    out = []
    for i in np.arange(inc, limit, inc):
        for j in range(repeat):
            p = -np.log(ks_2samp(data1[feature].sample(n=i).tolist(), data2[feature].sample(n=i).tolist())[1])
            out.append(p)
    return out


def get_x(inc: int, limit: int, repeat: int, offset: float):
    """
    Create pair-wise x value for -ln(p) plot

    :param inc: int, increment (generally 5)
    :param limit: int, upper limit of the plot
    :param repeat: int, how many runs to be calculated per condition (generally 50)
    :param offset: float, offset applied to avoid different datasets overlapping
    :return: out: list, list of x values
    """
    out = []
    for i in np.arange(inc, limit, inc):
        for j in range(repeat):
            x = i+offset
            out.append(x)
    return out


def plot_volcano_hit(pd_p: pd.DataFrame, pd_value: pd.DataFrame, pd_gamma: pd.DataFrame, feature: str, center: float,
                     threshold: float, show_gene: str, save_path: str):
    """
    Plot volcano plot for FRAP screen (with hit gene labeled)

    :param pd_p: pd.DataFrame, -ln(p) value dataframe
    :param pd_value: pd.DataFrame, average value dataframe
    :param pd_gamma: pd.DataFrame, score gamma (value-center)*[-ln(p)] dataframe
    :param feature: str, plotting feature, column name
    :param center: float, center used for calculating gamma, get from volcano plot-sum
    :param threshold: float, threshold used to call hit
                      fold of value difference from center/ std(control values) * 3 (roughly -ln(0.05))
    :param show_gene: str, only accepts 'N' or 'Y', whether display gene name on volcano plot
    :param save_path: str, saving path
    :return:
    """
    plt.figure(figsize=(9, 6))

    ctrl_std = np.std(pd_value[pd_value['sample'].isin(['WT', 'NT', 'NTC1', 'NTC2', 'E3', 'E4'])][feature])
    ishit = pd_gamma[feature] / ctrl_std >= threshold
    isgene = ~pd_gamma['sample'].isin(['WT', 'NT', 'NTC1', 'NTC2', 'E3', 'E4'])
    isNT = pd_gamma['sample'].isin(['NT', 'NTC1', 'NTC2', 'E3', 'E4'])
    isWT = pd_gamma['sample'].isin(['WT'])

    plt.scatter(pd_value[~ishit & isgene][feature], pd_p[~ishit & isgene][feature], s=6, c='#A9A9A9',
                label='gene non-hit')
    plt.scatter(pd_value[ishit & isgene][feature], pd_p[ishit & isgene][feature], s=6, c='#DC143C', label='gene hit')
    plt.scatter(pd_value[~ishit & isNT][feature], pd_p[~ishit & isNT][feature], s=6, c='#FFD700', label='NT non-hit')
    plt.scatter(pd_value[ishit & isNT][feature], pd_p[ishit & isNT][feature], s=6, c='#FF8C00', label='NT hit')
    # plt.scatter(pd_value[~ishit & isWT][feature], pd_p[~ishit & isWT][feature], s=6, c='#DDA0DD', label='WT non-hit')
    # plt.scatter(pd_value[ishit & isWT][feature], pd_p[ishit & isWT][feature], s=6, c='#8A2BE2', label='WT hit')

    ymax = np.ceil(max(pd_p[feature])) * 1.02
    xmin = min(pd_value[feature]) * 0.95
    xmax = max(pd_value[feature]) * 1.1

    plt.plot(np.linspace(xmin, xmax, 1000),
             np.abs(threshold / np.linspace((xmin-center) / ctrl_std, (xmax-center) / ctrl_std, 1000)), 'k--', lw=.5)

    if show_gene == 'Y':
        x_lst = pd_value[ishit & isgene][feature].tolist()
        y_lst = pd_p[ishit & isgene][feature].tolist()
        gene_lst = pd_value[ishit & isgene]['sample'].tolist()
        for i in range(len(x_lst)):
            plt.text(x_lst[i], y_lst[i], gene_lst[i], fontsize=6)

    plt.xlim((xmin, xmax))
    plt.ylim((0, ymax))
    plt.xlabel(feature)
    plt.ylabel('-ln(p)')
    plt.legend(loc=4, bbox_to_anchor=(0.7, 0, 0.3, 0.3))

    plt.savefig('%s%s_hit.pdf' % (save_path, feature))
    plt.close()