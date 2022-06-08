import numpy as np
from skimage.filters import threshold_otsu, threshold_local
from skimage.morphology import binary_dilation, binary_erosion, disk
from skimage.segmentation import clear_border
import shared.objects as obj
from scipy import ndimage
import math
import shared.dataframe as dat
from scipy.optimize import curve_fit

"""
# ---------------------------------------------------------------------------------------------------
# FUNCTIONS related with MATH
# ---------------------------------------------------------------------------------------------------

auto_correlation
    FUNCTION: auto_correlation function rewritten from matlab code contributed by Sarah Veatch
    SYNTAX:   auto_correlation(img: np.array, mask: np.array, rmax: int)

cart2pol
    FUNCTION: transform Cartesian coordinate into polar coordinates
    SYNTAX:   cart2pol(x: float, y: float)

cart2pol_matrix
    FUNCTION: Transform Cartesian coordinate into polar coordinates, matrix format
    SYNTAX:   cart2pol_matrix(x: np.array, y: np.array)

linear
    EQUATION: y = a * x + b
    SYNTAX:   linear(x, a, b)
    
r_square
    FUNCTION: calculate r_squared for a given curve fitting
    SYNTAX:   r_square(y: list, y_fit: list)

fitting_linear
    FUNCTION: perform linear fitting
    SYNTAX:   fitting_linear(x: list, y: list)  
"""


def auto_correlation(img: np.array, mask: np.array, rmax: int):
    """
    Auto_correlation function rewritten from matlab code contributed by Sarah Veatch

    :param img: np.array, color image for auto correlation function calculation
    :param mask: np.array, 0-1 image, segmented mask image
    :param rmax: int, padding range
    :return:
    """

    img_int = np.sum(img*mask)
    mask_int = np.sum(mask)

    (l1, l2) = img.shape

    mask_fft = np.real(np.fft.fftshift(np.real(np.fft.ifft2(np.abs(np.fft.fft2(dat.matrix_pad_with_zeros(mask, rmax, rmax))) ** 2))))
    img_fft = np.real(np.fft.fftshift(np.real(np.fft.ifft2(np.abs(np.fft.fft2(dat.matrix_pad_with_zeros(img * mask, rmax, rmax))) ** 2))))
    mask_fft[mask_fft == 0] = 10 ** -16
    img_fft[img_fft == 0] = 10 ** -16

    g1 = (mask_int ** 2)/(img_int ** 2)*img_fft/mask_fft
    gg = g1[math.floor((l1+rmax)/2+1)-rmax-1: math.floor((l1+rmax)/2+1)+rmax, math.floor((l2+rmax)/2+1)-rmax-1: math.floor((l2+rmax)/2+1)+rmax]

    xvals = np.array([np.ones(2*rmax+1)]).transpose() * np.arange(-rmax, rmax+1, 1)
    yvals = np.array([np.arange(-rmax, rmax+1, 1)]).transpose() * np.ones(2*rmax+1)
    zvals = gg

    r, theta = cart2pol_matrix(xvals, yvals)
    v = zvals

    ar = r.reshape(1, (2*rmax+1)**2)
    avals = v.reshape(1, (2*rmax+1)**2)

    ind = np.argsort(ar)
    rr = np.sort(ar)

    vv = np.array([avals[0][i] for i in ind])

    r = np.arange(0, math.ceil(max(rr[0]))+1, 1)
    n, _ = np.histogram(rr, bins=[i-0.5 for i in r]+[r[-1]+0.5])
    bin = np.digitize(rr, bins=[i-0.5 for i in r]+[r[-1]+0.5])

    g = []
    dg = []
    for j in range(rmax+1):
        m = (bin == j+1)*bin
        m = np.array([[i/(j+1) for i in m[0]]])
        n2 = sum(m[0])
        if n2 != 0:
            temp = m*vv
            g_temp = sum(temp[0])/n2
            g.append(g_temp)

            temp = m*np.array([[i-g_temp for i in vv[0]]]) ** 2
            dg_temp = math.sqrt(sum(temp[0]))/n2
            dg.append(dg_temp)

    r = np.arange(0, rmax+1, 1)
    gg[rmax+1][rmax+1] = 0

    return gg, r, g, dg


def cart2pol(x: float, y: float):
    """
     Transform Cartesian coordinate into polar coordinates

    :param x: float, cartesian coordinate x
    :param y: float, cartesian coordinate y
    :return:
    """
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)

    return r, theta


def cart2pol_matrix(x: np.array, y: np.array):
    """
    Transform Cartesian coordinate into polar coordinates, matrix format

    :param x: np.array, matrix of cartesian coordinates x
    :param y: np.array, matrix of cartesian coordinates y
    :return:
    """
    r = np.zeros_like(x, dtype=float)
    theta = np.zeros_like(x, dtype=float)
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            r[i][j], theta[i][j] = cart2pol(x[i][j], y[i][j])

    return r, theta


# linear regression
def linear(x, a, b):
    """
        Linear regression function
    """

    return a * x + b


def linear_b0(x, a):
    return a * x


# r_square calculation
def r_square(y: list, y_fit: list):
    """
    calculate r_squared for a given curve fitting

    :param y: values before fitting
    :param y_fit: values from fitting
    :return: r_squared
    """

    ss_res = np.sum([(a - b) ** 2 for a, b in zip(y, y_fit)])
    ss_tot = np.sum([(a - np.mean(y)) ** 2 for a in y])
    r2 = 1 - (ss_res / ss_tot)
    return r2


def fitting_linear(x: list, y: list):
    """
    Perform linear fitting

    y = a * x + b

    :param x:
    :param y:
    :return: y_fit:
             r2:
             a:
             b:
    """

    try:
        popt, _ = curve_fit(linear, x, y)
        a, b = popt
    except RuntimeError:
        a = np.nan
        b = np.nan

    y_fit = []
    for j in range(len(y)):
        y_fit.append(linear(x[j], a, b))
    r2 = r_square(y, y_fit)

    return y_fit, r2, a, b


def fitting_linear_b0(x: list, y: list):
    """
    Perform linear fitting

    y = a * x

    :param x:
    :param y:
    :return: y_fit:
             r2:
             a:
             b:
    """

    try:
        a, _ = curve_fit(linear_b0, x, y)
    except RuntimeError:
        a = np.nan

    y_fit = []
    for j in range(len(y)):
        y_fit.append(linear_b0(x[j], a))
    r2 = r_square(y, y_fit)

    return y_fit, r2, a
