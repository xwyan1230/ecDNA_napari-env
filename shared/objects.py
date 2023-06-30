import numpy as np
from skimage.morphology import remove_small_objects, medial_axis, extrema, binary_dilation, dilation
from skimage.measure import label, regionprops
from skimage.filters import sobel
from skimage import segmentation
import math
import random

"""
# ---------------------------------------------------------------------------------------------------
# FUNCTIONS for 0-AND-1 NP.ARRAY (BINARY IMAGE)
# ---------------------------------------------------------------------------------------------------

remove_small
    FUNCTION: remove objects smaller than the specified size
    SYNTAX:   remove_small(obj: np.array, min_size=10)

remove_large
    FUNCTION: remove objects larger than the specified size
    SYNTAX:   remove_large(obj: np.array, max_size=1000)

get_intensity
    FUNCTION: measure mean intensity time series for all given objects
    SYNTAX:   get_intensity(obj: np.array, pixels_tseries: list)

label_watershed
    FUNCTION: separate touching objects based on distance map (similar to imageJ watershed)
    SYNTAX:   label_watershed(obj: np.array)
    
label_remove_large
    FUNCTION: remove objects larger than the specified size from labeled image
    SYNTAX:   label_remove_large(label_obj: np.array, max_size: int)

label_remove_small
    FUNCTION: remove objects smaller than the specified size from labeled image
    SYNTAX:   label_remove_small(label_obj: np.array, min_size: int)
    
label_remove_low_circ
    FUNCTION: remove objects whose circularity are smaller than the specified threshold from labeled image
    SYNTAX:   label_remove_low_circ(label_obj: np.array, thresh: float)
    
label_watershed
    FUNCTION: separate touching objects based on distance map (similar to imageJ watershed)
    SYNTAX:   label_watershed(obj: np.array, maxima_threshold)

label_resort
    FUNCTION: re-label random labeled image into sequential labeled image
    SYNTAX:   label_resort(label_obj: np.array)

object_count
    FUNCTION: count the number of objects in given image
    SYNTAX:   object_count(obj: np.array)
"""


def remove_small(obj: np.array, min_size=10):
    """
    Remove objects smaller than the specified size.

    Expects obj to be an integer image array with objects labeled with 1, and removes objects
    smaller than min_size.

    :param obj: np.array, 0-and-1
    :param min_size: int, optional (default: 10)
                The smallest allowable object size.
    :return: out: nd.array, 0-and-1, same shape and type as input obj

    """

    obj_bool = np.array(obj, bool)
    obj_mask = remove_small_objects(obj_bool, min_size)
    out = np.zeros_like(obj)
    out[obj_mask] = 1

    return out


def remove_large(obj: np.array, max_size=1000):
    """
    Remove objects larger than the specified size.

    Expects obj to be an integer image array with objects labeled with 1, and removes objects
    larger than max_size.

    :param obj: np.array, 0-and-1
    :param max_size: int, optional (default: 1000)
                The largest allowable object size.
    :return: out: np.array, 0-and-1, same shape and type as input obj
    """

    obj_bool = np.array(obj, bool)
    obj_mask = remove_small_objects(obj_bool, max_size)
    out = obj.copy()
    out[obj_mask] = 0

    return out


def get_intensity(label_obj: np.array, pixels_tseries: list):
    """
    Measure mean intensity time series for all given objects

    Usage examples:
    1) measure bleach spots/ctrl spots intensity series

    :param label_obj: np.array, 0-and-1 object mask
    :param pixels_tseries: list, pixels time series
                e.g. [pixels_t0, pixels_t1, pixels_t2, ...]
    :return: obj_int_tseries: list
                list of intensity time series
    """

    max_t = len(pixels_tseries)
    obj_int_tseries = [[] for _ in range(len(np.unique(label_obj))-1)]

    for t in range(0, max_t):
        # measure mean intensity for objects
        obj_props = regionprops(label_obj, pixels_tseries[t])
        for i in range(len(obj_props)):
            obj_int_tseries[i].append(obj_props[i].mean_intensity)

    return obj_int_tseries


def label_watershed(obj: np.array, maxima_threshold):
    """
    Separate touching objects based on distance map (similar to imageJ watershed)

    :param obj: np.array, 0-and-1
    :param maxima_threshold: threshold for identify maxima
    :return: seg: np.array, grey scale with different objects labeled with different numbers
    """
    _, dis = medial_axis(obj, return_distance=True)
    maxima = extrema.h_maxima(dis, maxima_threshold)
    # maxima_threshold for Jess data = 1
    # maxima_threshold for Jose 60x data = 10
    # maxima_threshold for Jose 40x data = 20
    maxima_mask = binary_dilation(maxima)
    for i in range(6):
        maxima_mask = binary_dilation(maxima_mask)

    label_maxima = label(maxima_mask, connectivity=2)
    markers = label_maxima.copy()
    markers[obj == 0] = np.amax(label_maxima) + 1
    elevation_map = sobel(obj)
    label_obj = segmentation.watershed(elevation_map, markers)
    label_obj[label_obj == np.amax(label_maxima) + 1] = 0

    return label_obj


def label_remove_large(label_obj: np.array, max_size: int):
    """
    Remove objects larger than the specified size from labeled image

    :param label_obj: np.array, grey scale labeled image
    :param max_size: int
                The largest allowable object size
    :return: out: np.array, grey scale labeled image with large objects removed
    """
    out = np.zeros_like(label_obj)
    obj_prop = regionprops(label_obj)
    for i in range(len(obj_prop)):
        if obj_prop[i].area <= max_size:
            out[label_obj == obj_prop[i].label] = obj_prop[i].label

    return out


def label_remove_small(label_obj: np.array, min_size: int):
    """
    Remove objects smaller than the specified size from labeled image

    :param label_obj: np.array, grey scale labeled image
    :param min_size: int
                The smallest allowable object size
    :return: out: np.array, grey scale labeled image with small objects removed
    """
    out = np.zeros_like(label_obj)
    obj_prop = regionprops(label_obj)
    for i in range(len(obj_prop)):
        if obj_prop[i].area >= min_size:
            out[label_obj == obj_prop[i].label] = obj_prop[i].label

    return out


def label_remove_small_large_resort(label_obj: np.array, min_size: int, max_size: int):
    """
        Remove objects smaller or larger than the specified size from labeled image and sort objects

        :param label_obj: np.array, grey scale labeled image
        :param min_size: int
                    The smallest allowable object size
        :param max_size: int
                    The largest allowable object size
        :return: out: np.array, grey scale labeled image with small objects removed
        """
    count = 1
    out = np.zeros_like(label_obj)
    obj_prop = regionprops(label_obj)
    for i in range(len(obj_prop)):
        """if i % 100 == 0:
            print("%s/%s" % (i, len(obj_prop)))"""
        if (obj_prop[i].area >= min_size) & (obj_prop[i].area <= max_size):
            out[label_obj == obj_prop[i].label] = count
            count = count + 1

    return out


def label_remove_low_circ(label_obj: np.array, thresh: float):
    """
    Remove objects whose circularity are smaller than the specified threshold from labeled image

    Note: this step will also remove any objects whose size are smaller than 50

    :param label_obj: np.array, grey scale labeled image
    :param thresh: float
                The smallest allowable circularity
    :return: out: np.array, grey scale labeled image with small circularity objects removed
    """
    obj_prop = regionprops(label_obj)
    out = np.zeros_like(label_obj, dtype=float)
    n = 1
    for i in obj_prop:
        if i.area >= 50:
            circ = (4 * math.pi * i.area) / (i.perimeter ** 2)
            if circ > thresh:
                out[label_obj == i.label] = n
                n = n+1
    return out


def label_watershed(obj: np.array, maxima_threshold):
    """
    Separate touching objects based on distance map (similar to imageJ watershed)

    :param obj: np.array, 0-and-1
    :param maxima_threshold: threshold for identify maxima
    :return: seg: np.array, grey scale with different objects labeled with different numbers
    """
    _, dis = medial_axis(obj, return_distance=True)
    maxima = extrema.h_maxima(dis, maxima_threshold)
    # maxima_threshold for Jess data = 1
    # maxima_threshold for Jose 60x data = 10
    # maxima_threshold for Jose 40x data = 20
    maxima_mask = binary_dilation(maxima)
    for i in range(6):
        maxima_mask = binary_dilation(maxima_mask)

    label_maxima = label(maxima_mask, connectivity=2)
    markers = label_maxima.copy()
    markers[obj == 0] = np.amax(label_maxima) + 1
    elevation_map = sobel(obj)
    label_obj = segmentation.watershed(elevation_map, markers)
    label_obj[label_obj == np.amax(label_maxima) + 1] = 0

    return label_obj


def label_resort(label_obj: np.array):
    """
    Re-label random labeled image into sequential labeled image.

    :param label_obj: np.array, grey scale labeled image
    :return: label_out: np.array, sorted grey scale labeled image
    """
    count = 1
    label_out = np.zeros_like(label_obj)
    for i in range(np.amax(label_obj)):
        temp = np.zeros_like(label_obj)
        temp[label_obj == i + 1] = 1
        if object_count(temp) > 0:
            label_out[label_obj == i + 1] = count
            count = count + 1

    return label_out


def label_resort_print(label_obj: np.array):
    """
    Re-label random labeled image into sequential labeled image.

    :param label_obj: np.array, grey scale labeled image
    :return: label_out: np.array, sorted grey scale labeled image
    """
    print("label resort ongoing...")
    count = 1
    label_out = np.zeros_like(label_obj)
    for i in range(np.amax(label_obj)):
        temp = np.zeros_like(label_obj)
        temp[label_obj == i + 1] = 1
        if object_count(temp) > 0:
            label_out[label_obj == i + 1] = count
            count = count + 1
    print("label resort done!")
    return label_out


def object_count(obj: np.array):
    """
    Count the number of objects in given image.

    :param obj: np.array, 0-and-1
    :return: count_obj: number of objects.
    """
    label_obj = label(obj, connectivity=1)
    obj_prop = regionprops(label_obj)
    count_obj = len(obj_prop)

    return count_obj


