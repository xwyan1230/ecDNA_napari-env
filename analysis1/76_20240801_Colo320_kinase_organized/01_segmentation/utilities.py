import numpy as np
from skimage.measure import regionprops, label
from skimage.filters import threshold_otsu
from skimage.morphology import medial_axis, extrema, binary_dilation, dilation, erosion, disk
from skimage.filters import sobel
from skimage.segmentation import clear_border, watershed
import math


def img_crop(img, i):
    if i == 0:
        out = img[0:1440, 0:1920]
    elif (i == 1) | (i == 2):
        out = img[0:1440, 576:1920]
    elif (i == 3) | (i == 6):
        out = img[432:1440, 0:1920]
    elif (i == 4) | (i == 5):
        out = img[432:1440, 0:1344]
    elif (i == 7) | (i == 8):
        out = img[432:1440, 576:1920]
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
        if (obj_prop[i].area >= min_size) & (obj_prop[i].area <= max_size):
            out[label_obj == obj_prop[i].label] = count
            count = count + 1
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
    maxima_mask = binary_dilation(maxima)
    for i in range(6):
        maxima_mask = binary_dilation(maxima_mask)

    label_maxima = label(maxima_mask, connectivity=2)
    markers = label_maxima.copy()
    markers[obj == 0] = np.amax(label_maxima) + 1
    elevation_map = sobel(obj)
    label_obj = watershed(elevation_map, markers)
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


def cell_seg_fluorescent(img: np.array, otsu_factor=1.5, maxima_threshold=1, max_size=1800,
                         min_size=300, circ_thresh=0.6, threshold=0):
    """
    Perform cell segmentation from a fluorescent image

    Algorithm description:
    Otsu segmentation followed by watershed and filtering

    :param img: np.array
                    fluorescent image
    :param otsu_factor: float
                    factor multiplied to otsu threshold used to perform otsu thresholding
    :param maxima_threshold: int
                    threshold for identify maxima during watershed
    :param max_size: float
                    maximum allowable object size
    :param min_size: float
                    minimum allowable object size
    :param circ_thresh: float
                    minimum allowable circularity
    :param threshold: artificial threshold for segmentation
    :return:
    """
    if threshold == 0:
        thresh_val = threshold_otsu(img)
        out = img > thresh_val * otsu_factor
    else:
        out = img > threshold
    out = label_watershed(out, maxima_threshold)
    out = clear_border(out)
    out = label_remove_large(out, max_size)
    out = label_remove_small(out, min_size)
    out = label_remove_low_circ(out, circ_thresh)
    out = dilation(out, disk(3))
    out = erosion(out, disk(3))

    return out


def img_factor(i):
    if i == 0:
        out = 1
    elif (i == 1) | (i == 2) | (i == 3) | (i == 6):
        out = 0.7
    elif (i == 4) | (i == 5) | (i == 7) | (i == 8):
        out = 0.49
    return out