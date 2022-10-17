import numpy as np
from skimage.filters import threshold_otsu, threshold_local, threshold_yen, sobel
from skimage.segmentation import clear_border
import shared.objects as obj
from scipy import ndimage
from skimage import segmentation
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk
from shared.objects import remove_large, remove_small
from skimage.measure import label, regionprops_table, regionprops
import pandas as pd
import shared.image as ima
import napari


"""
# ---------------------------------------------------------------------------------------------------
# FUNCTIONS for OBJECT SEGMENTATION
# ---------------------------------------------------------------------------------------------------

cell_seg_trans
    FUNCTION: perform cell segmentation from a transmitted light image
    SYNTAX:   cell_seg_trans(img_trans: np.array, max_size=1000, min_size=100, selem_num=8)
    
cell_seg_fluorescent   
    FUNCTION: perform cell segmentation from a fluorescent image
    SYNTAX:   cell_seg_fluorescent(img: np.array, otsu_factor=1.5, maxima_threshold=1, max_size=1800, 
                                   min_size=300, circ_thresh=0.6)

nuclear_seg
    FUNCTION: perform nuclear segmentation from a fluorescent image
    SYNTAX:   nuclear_seg(img: np.array, local_factor=99, clearance_threshold=300, maxima_threshold=10, min_size=4000, 
              max_size=25000)

segment_watershed
    FUNCTION: returns an np.array that is the segmented version of the input
    SYNTAX:   segment_watershed(pixels: np.array, extreme_val: int, bg_val: int)

find_blobs
    FUNCTION: find blobs in image
    SYNTAX:   find_blobs(pixels: np.array, binary_global: np.array, extreme_val: int, bg_val=200, 
              max_size=1000)

get_binary_global
    FUNCTION: calculate binary global thresholding image
    SYNTAX:   get_binary_global(pixels: np.array, threshold_method='na', min_size=5, max_size=1000)

find_organelle
    FUNCTION: find organelle (nucleoli or SG) from a given image
    SYNTAX:   find_organelle(pixels: np.array, global_thresholding='na', extreme_val=500, bg_val=200,
              min_size=5, max_size=1000)

get_bg_int
    FUNCTION: measure background intensities from a given movie
    SYNTAX:   get_bg_int(pixels_tseries: list)

get_mean_int_of_max_area
    FUNCTION: get mean intensity of maximum area object
    SYNTAX:   get_mean_int_of_max_area(props: regionprops)

obj_to_convex
    FUNCTION: transform objects into its corresponding convex objects
    SYNTAX:   obj_to_convex(pixels: np.array)

obj_to_convex_filter
    FUNCTION: transform objects into its corresponding convex objects, only real area/convex area larger 
              than threshold will be converted
    SYNTAX:   obj_to_convex_filter(img_obj: np.array, threshold=0.9)
    
filter_solidity
    FUNCTION: filter labeled objects based on solidity
    SYNTAX:   filter_solidity(pixels: np.array, threshold=0.9)

filter_mean_int
    FUNCTION: filter labeled objects based on mean sub_obj intensity
    SYNTAX:   filter_mean_int(img_obj: np.array, img_sub_obj: np.array, img: np.array, threshold: float)

get_bg_img
    FUNCTION: generate background image using otsu
    SYNTAX:   get_bg_img(img: np.array)
"""


def cell_seg_trans(img_trans: np.array, max_size=1000, min_size=100, selem_num=8):
    """
    Perform cell segmentation from a transmitted light image

    Algorithm description:
    Segment cells by identifying highlighted parts surrounded by dark boundaries. This algorithm will
    potentially miss lots of cells with possible false-positives of spaces between cells which is
    around similar size of a given cell.
    Perform otsu thresholding to identify high intensity region, serial erosion and dilation to better
    close cell regions, size selection to exclude large background region and small debris.
    This method has only been tested for images acquired with EVOS-M5000.

    :param img_trans: np.array
                        transmitted light image
    :param max_size: int, optional, default = 1000
                        maximum size of a cell
    :param min_size: int, optional, default = 100
                        minimum size of a cell
    :param selem_num: int, optional, default = 8
                        number of disk size, tested from 5-9 and 7-9 seem all fine
    :return: out: np.array, binary image
                same size as original transmitted light image
                1: cell region
                0: background
    """

    thresh_val = threshold_otsu(img_trans)
    out = img_trans > thresh_val

    selem = disk(selem_num)
    out = binary_erosion(out, selem)
    out = binary_dilation(out, selem)

    out = obj.remove_large(out, max_size=max_size)
    out = obj.remove_small(out, min_size=min_size)

    return out


def cell_seg_fluorescent(img: np.array, otsu_factor=1.5, maxima_threshold=1, max_size=1800,
                         min_size=300, circ_thresh=0.6):
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
    :return:
    """
    thresh_val = threshold_otsu(img)
    out = img > thresh_val * otsu_factor
    out = obj.label_watershed(out, maxima_threshold)
    out = obj.label_remove_large(out, max_size)
    out = obj.label_remove_small(out, min_size)
    out = obj.label_remove_low_circ(out, circ_thresh)

    return out


def nuclear_seg(img: np.array, local_factor=99, clearance_threshold=300, maxima_threshold=10, min_size=4000, max_size=25000):
    """
    Perform nuclear segmentation from a fluorescent image

    tested by Paul Mischel Leica Scope

    :param img: np.array
                    fluorescent image
    :param local_factor: int, odd number
                    factor used to perform local thresholding
                    for ColoDM under Paul Mischel Leica scope, 99
    :param clearance_threshold: int
                    threshold used to clear background
                    default: 300
    :param maxima_threshold: int
                    threshold used in label_watershed
                    for ColoDM under Paul Mischel Leica scope, 10
    :param min_size: int
                    minimum allowable object size
    ;param max_size: int
                    maximum allowable object size
    :return: out: np.array
                    labeled nuclear img
    """
    # global thresholding to determine rough location of nuclei
    global_threshold_val = threshold_otsu(img)
    # determine background region
    bg = img > global_threshold_val
    # perform local thresholding to identify nuclei
    local = threshold_local(img, local_factor)
    out = img > local
    # clear background
    out[bg == 0] = 0
    # eliminate nuclei that touching boundary
    out = clear_border(out)
    # one round of erosion/dilation and clearance to clear background
    out = binary_erosion(out)
    out = obj.remove_small(out, clearance_threshold)
    out = binary_dilation(out)
    # fill nuclei holes
    out = ndimage.binary_fill_holes(out)
    # separate touching nuclei
    out = obj.label_watershed(out, maxima_threshold)
    # filter smaller objects
    out = obj.label_remove_small(out, min_size)
    out = obj.label_remove_large(out, max_size)
    out = obj.label_resort(out)

    return out


def nuclear_seg1(img: np.array, local_factor=99, clearance_threshold=300, maxima_threshold=10, min_size=4000, max_size=25000):
    """
    Perform nuclear segmentation from a fluorescent image

    tested by Paul Mischel Leica Scope

    :param img: np.array
                    fluorescent image
    :param local_factor: int, odd number
                    factor used to perform local thresholding
                    for ColoDM under Paul Mischel Leica scope, 99
    :param clearance_threshold: int
                    threshold used to clear background
                    default: 300
    :param maxima_threshold: int
                    threshold used in label_watershed
                    for ColoDM under Paul Mischel Leica scope, 10
    :param min_size: int
                    minimum allowable object size
    ;param max_size: int
                    maximum allowable object size
    :return: out: np.array
                    labeled nuclear img
    """
    # global thresholding to determine rough location of nuclei
    global_threshold_val = threshold_otsu(img)
    # determine background region
    bg = img > global_threshold_val
    # perform local thresholding to identify nuclei
    local = threshold_local(img, local_factor)
    out = img > local
    # clear background
    out[bg == 0] = 0
    # one round of erosion/dilation and clearance to clear background
    out = binary_erosion(out)
    out = binary_dilation(out, disk(4))
    out = obj.remove_small(out, clearance_threshold)
    # fill nuclei holes
    out = ndimage.binary_fill_holes(out)
    # eliminate nuclei that touching boundary
    out = clear_border(out)
    out = binary_erosion(out, disk(3))
    # separate touching nuclei
    out = obj.label_watershed(out, maxima_threshold)
    # filter smaller objects
    out = obj.label_remove_small(out, min_size)
    out = obj.label_remove_large(out, max_size)
    out = obj.label_resort(out)

    return out


def segment_watershed(pixels: np.array, extreme_val: int, bg_val: int):
    """
    Returns an np.array that is the segmented version of the input

        Finds local maxima using extreme_val as the minimal height for the maximum,
        Then uses the local maxima and background pixels (pixels that are
        smaller than bg_val) to execute a watershed segmentation
    """
    maxima = extrema.h_maxima(pixels, extreme_val)
    elevation_map = sobel(pixels)
    markers = np.zeros_like(pixels)
    markers[pixels < bg_val] = 1
    markers[maxima == 1] = 2

    return segmentation.watershed(elevation_map, markers)


def find_blobs(pixels: np.array, binary_global: np.array, extreme_val: int, bg_val=200, max_size=1000):
    """
    Find "blobs" in image.  Current strategy: use a segmentation.watershed on the input pixes,
    using extreme_val to find local maxima, adn bg_val to find background.  Combine the
    watershed with a globally threshold image (using logical OR) binary_global.

    :param pixels: input image
    :param binary_global: binary threshold image gain from global thresholding
    :param extreme_val: used to find local maxima
    :param bg_val: used to define background for watershed
    :param max_size: maximum size of the blobs
    :return: segmented image of same size as input

    """
    if np.amax(pixels) < 1000:
        merge = np.zeros_like(pixels)
    else:
        seg_wat = segment_watershed(pixels, extreme_val, bg_val)
        merge = np.zeros_like(pixels)
        merge[seg_wat == 2] = 1
        merge = remove_large(merge, max_size)
        binary_global = remove_large(binary_global, max_size)
        merge[binary_global == 1] = 1

    return merge


def get_binary_global(pixels: np.array, threshold_method='na', min_size=5, max_size=1000, local_param=21):
    """
    Calculate binary global thresholding image

    :param pixels: np.array
    :param threshold_method: method used to perform global thresholding, enable 'na',
                'otsu', 'yen', 'local-nucleoli' and 'local-sg', 'local-sg1'
                'na': not applied, return a black image
                'otsu': otsu thresholding + one round of erosion/dilation
                'yen': yen thresholding + one round of erosion/dilation
                'local-nucleoli': otsu & local thresholding for nucleoli identification
                'local-sg': otsu & local thresholding for stress granule identification
    :param min_size: minimum size of blobs
    :param max_size: maximum size of blobs
    :param local_param: parameter for local thresholding
    :return: out: 0-and-1 np.array, binary global thresholding image

    """

    check_lst = ['na', 'otsu', 'yen', 'local-nucleoli', 'local-sg', 'local-sg1']
    if threshold_method not in check_lst:
        raise ValueError("global thresholding method only accepts %s. Got %s" % (check_lst, threshold_method))

    elif (threshold_method == 'otsu') | (threshold_method == 'yen'):
        if threshold_method == 'otsu':
            global_threshold_val = threshold_otsu(pixels)
            # Threshold value to create global threshold.  try: threshold_otsu(pixels)
            # 0: does not apply global thresholding
        else:
            global_threshold_val = threshold_yen(pixels)
        out = pixels > global_threshold_val
        # one round of erosion/dilation to clear out boundary
        out = binary_erosion(out)
        out = binary_dilation(out)

    elif threshold_method == 'local-nucleoli':
        # use otsu thresholding to determine background region
        global_threshold_val = threshold_otsu(pixels)
        bg = pixels > global_threshold_val
        # apply local thresholding
        local = threshold_local(pixels, local_param)  # 21: specific for nucleoli
        out = pixels > local
        # remove large connected areas
        out = obj.remove_large(out, max_size)
        # combine with otsu thresholding to determine background region
        out[bg == 0] = 0
        # two rounds of erosion/dilation and remove_small to clear out background
        out = binary_erosion(out)
        out = obj.remove_small(out, min_size)
        out = binary_dilation(out)

    elif threshold_method == 'local-sg':
        # use otsu thresholding to determine background region
        global_threshold_val = threshold_otsu(pixels)
        bg = pixels > global_threshold_val
        # apply local thresholding
        local = threshold_local(pixels, local_param)  # 21: specific for nucleoli
        out = pixels > local
        # remove large connected areas
        out = obj.remove_large(out, max_size)
        # combine with otsu thresholding to determine background region
        out[bg == 0] = 0
        out = ndimage.binary_fill_holes(out)

    elif threshold_method == 'local-sg1':
        # use otsu thresholding to determine background region
        global_threshold_val = threshold_otsu(pixels)
        global_threshold_val1 = threshold_yen(pixels)
        bg = pixels > global_threshold_val
        bg1 = pixels > global_threshold_val1
        bg2 = pixels > 2500
        # apply local thresholding
        local = threshold_local(pixels, 51)
        out = pixels > local
        # remove large connected areas
        out = obj.remove_large(out, max_size)
        # combine with otsu thresholding to determine background region
        out[bg == 0] = 0
        out[bg1 == 0] = 0
        out[bg2 == 0] = 0

    else:
        out = np.zeros_like(pixels)

    return out


def find_organelle(pixels: np.array, global_thresholding='na', extreme_val=500, bg_val=200,
                   min_size=5, max_size=1000, local_param=21):
    """
    Find organelle (nucleoli or SG) from a given image.

    Expects pixels to be an array, and finds nucleoli objects using watershed by
    flooding approach with indicated global thresholding methods (supports 'na',
    'otsu', 'yen', 'local-nucleoli' and 'local-sg', 'local-sg1').  Founded organelles
    are filtered by default location filter (filter out organelles located at the
    boundary of the image) and size filter.

    :param pixels: np.array (non-negative int type)
                Image pixel
    :param global_thresholding: only accepts 'na', 'otsu', 'yen', 'local-nucleoli'
                or 'local-sg'
                optional (default: 'na')
                Whether or not ('na') to apply global thresholding method and
                which method to apply
                nucleoli default: 'local-nucleoli'
                SG default: 'local-sg'
    :param extreme_val: int, optional (default: 500)
                Used in shared.find_blobs.segment_watershed to find local maxima
    :param bg_val: int, optional (default: 200)
                Used in shared.find_blobs.segment_watershed to define background
                for watershed
    :param min_size: int, optional (default: 5)
                The smallest allowable organelle size.
                nucleoli default: 10
                SG default: 5
    :param max_size: int, optional (default: 1000)
                The largest allowable organelle size.
                nucleoli default: 1000
                SG default: 350
    :param local_param: int, optional (default: 21)
                parameter for local thresholding
    :returns nucleoli_filtered: 0-and-1 np.array, same shape and type as input img
                Binary array with found nucleoli labeled with 1.
    """

    # Check global thresholding options
    # Raise value error if not 'na', 'otsu' or 'yen'
    check_lst = ['na', 'otsu', 'yen', 'local-nucleoli', 'local-sg', 'local-sg1']
    if global_thresholding not in check_lst:
        raise ValueError("global thresholding method only accepts %s. Got %s" % (check_lst, global_thresholding))

    # find organelle
    organelle = find_blobs(pixels, get_binary_global(pixels, global_thresholding, min_size, max_size, local_param),
                           extreme_val, bg_val, max_size)

    # Filters:
    # Location filter: remove artifacts connected to image border
    organelle_filtered = clear_border(organelle)

    # Size filter: default [10,1000]
    organelle_filtered = remove_small(organelle_filtered, min_size)
    organelle_filtered = remove_large(organelle_filtered, max_size)

    return organelle, organelle_filtered


def get_bg_int(pixels_tseries: list):
    """
    Measure background intensities from a given movie

    Algorithm description:
    For each time frame from the movie pixels_tseries, generate binary image of regions from otsu
    thresholding, remove regions whose area < 50 and return the mean intensity of the largest area
    as background intensity for this given time frame.

    Usage examples:
    1) used for background correction

    Note:
    1) function was originally designed for getting background series from a movie, but can also
        be applied to get background intensity from a single image. If so, please do:
        bg_int = get_bg_int([pixels])[0]

    :param pixels_tseries: list
                time series of np.array (movie), e.g. [pixel1, pixel2, ...]
                pixels_tseries[i]: np.array, pixels at time frame i
    :return: bg_int_tseries: list, 0 indicates background intensity detection failure
                list of background intensities, e.g. [bg_1, bg_2, ...]
                t_bg_int[i]: float, bg_int at frame i

    """
    bg_int_tseries = []
    for i in range(len(pixels_tseries)):
        # get regions based on otsu thresholding
        bg = pixels_tseries[i] < threshold_otsu(pixels_tseries[i])
        # smooth region
        bg = binary_dilation(bg)
        bg = binary_dilation(bg)
        bg = binary_erosion(bg)
        bg = binary_erosion(bg)
        # remove regions < 50
        bg = obj.remove_small(bg, 50)
        # measure bg object properties
        bg_props = regionprops_table(label(bg), pixels_tseries[i], properties=['area', 'mean_intensity'])
        bg_prop_pd = pd.DataFrame(bg_props)

        # high intensity image, do not have pixel intensity of any connected 50 pixels < bg_thresh
        if len(bg_prop_pd) == 0:
            # set bg_int as 0 to get the script going without interruption
            # 0 should not affect bg intensity curve fitting
            bg_int_tseries.append(0)
        elif len(bg_prop_pd) == 1:
            bg_int_tseries.append(bg_prop_pd['mean_intensity'][0])
        else:
            # find the mean_intensity of the largest area
            max_area = 0
            bg_int_temp = 0
            for j in range(len(bg_prop_pd)):
                if bg_prop_pd['area'][j] > max_area:
                    max_area = bg_prop_pd['area'][j]
                    bg_int_temp = bg_prop_pd['mean_intensity'][j]
            bg_int_tseries.append(bg_int_temp)

    return bg_int_tseries


def get_mean_int_of_max_area(props):
    """
    Get mean intensity of maximum area object

    :param props: regionprops
    :return:
    """
    max_area = 0
    mean_int = 0
    for i in range(len(props)):
        if props[i].area > max_area:
            max_area = props[i].area
            mean_int = props[i].mean_intensity

    return mean_int, max_area


def obj_to_convex(img_obj: np.array):
    """
    Transform objects into its corresponding convex objects
    :param img_obj: np.array, labeled image
    :return:
    """
    out = np.zeros_like(img_obj)
    props = regionprops(img_obj)
    for i in range(len(props)):
        convex_local = props[i].convex_image
        centroid = props[i].centroid
        centroid_convex = regionprops(label(convex_local))[0].centroid
        out = ima.image_paste_fix_value(out, convex_local, [int(centroid[0] - centroid_convex[0]),
                                                            int(centroid[1] - centroid_convex[1])], i)

    return out


def obj_to_convex_filter(img_obj: np.array, threshold=0.9):
    """
    Transform objects into its corresponding convex objects, only real area/convex area larger than threshold
    will be converted
    :param img_obj: np.array, labeled image
    :param threshold: used to filter ratio between real area/convex area, smaller than threshold will not be
                    converted into convex
    :return:
    """
    out = np.zeros_like(img_obj)
    props = regionprops(img_obj)
    for i in range(len(props)):
        convex_local = props[i].convex_image
        area_ratio = props[i].area/convex_local.sum()
        if area_ratio > threshold:
            centroid = props[i].centroid
            centroid_convex = regionprops(label(convex_local))[0].centroid
            out = ima.image_paste_fix_value(out, convex_local, [int(centroid[0] - centroid_convex[0]),
                                                                int(centroid[1] - centroid_convex[1])], i)

    return out


def filter_solidity(img_obj: np.array, threshold=0.9):
    """
    Filter labeled objects based on solidity
    :param img_obj: np.array, labeled image
    :param threshold: float, smaller than this number will be excluded
    :return:
    """
    out = np.zeros_like(img_obj)
    props = regionprops(img_obj)
    j = 1
    for i in range(len(props)):
        if props[i].solidity >= threshold:
            out[img_obj == props[i].label] = j
            j = j+1

    return out


def filter_mean_int(img_obj: np.array, img_sub_obj: np.array, img: np.array, threshold: float):
    """
    Filter labeled objects based on mean sub_obj intensity
    :param img_obj: np.array, labeled image, first tier (for example: nuclear segmentation)
    :param img_sub_obj: np.array, labeled image, second tier
            (for example: sub-nuclear segmentation, like ecDNA segmentation)
    :param img: np.array, fluorescent image for measuring intensity
    :param threshold: float, smaller than this number will be excluded
    :return:
    """
    out = np.zeros_like(img_obj)
    props = regionprops(img_obj)
    j = 1
    for i in range(len(props)):
        temp = img_sub_obj.copy()
        temp[img_obj != props[i].label] = 0
        props_temp = regionprops(label(temp), img)
        sum_intensity = 0
        sum_area = 0
        sub_obj_mean_intensity = 0
        for m in range(len(props_temp)):
            sum_intensity = sum_intensity + props_temp[m].mean_intensity*props_temp[m].area
            sum_area = sum_area + props_temp[m].area
            sub_obj_mean_intensity = sum_intensity*1.0/sum_area if sum_area != 0 else 0
        if sub_obj_mean_intensity >= threshold:
            out[img_obj == props[i].label] = j
            j = j+1

    return out


def get_bg_img(img: np.array):
    """
    Generate background image using otsu

    :param img: np.array, input image
    :return:
    """
    bg = np.ones_like(img)
    sample = img > threshold_otsu(img)
    sample = binary_dilation(sample)
    sample = binary_erosion(sample, disk(2))
    sample = binary_dilation(sample)
    bg[sample == 1] = 0

    return bg, sample


def puncta_seg(img: np.array, nuclear_seg: np.array, local_size: int, local_cycle_factor=15, filter='T',
               min_threshold_factor=6, min_size=4, min_threshold_first_round_factor=0.5, training='F'):
    """
    Segmentation for ecDNA or other puncta signal

    :param img: np.array, img to be segmented
    :param nuclear_seg: np.array, nuclear segmentation image
    :param local_size: int, for local image processing
    :param min_threshold_factor: int, minimum intensity threshold for puncta identification
    :param filter: str, only accepts 'T' or 'F'
    :param local_cycle_factor: int, cycle number for local segmentation
    :return:
    """
    nuclear_props = regionprops(nuclear_seg)
    img_seg = np.zeros_like(img)

    for i in range(len(nuclear_props)):
        original_centroid_nuclear = nuclear_props[i].centroid
        position = ima.img_local_position(nuclear_seg, original_centroid_nuclear, local_size)
        local_nuclear_seg_convex = ima.img_local_seg(nuclear_seg, position, nuclear_props[i].label)
        local_img = img.copy()
        local_img = local_img[position[0]:position[1], position[2]:position[3]]
        local_nuclear_props = regionprops(label(local_nuclear_seg_convex))
        local_centroid = local_nuclear_props[0].centroid

        # ecDNA segmentation
        local_img_singlet = local_img.copy()
        local_img_singlet[local_nuclear_seg_convex == 0] = 0
        otsu_threshold_val_local_img = threshold_otsu(local_img_singlet)

        if otsu_threshold_val_local_img == 0:
            print("skip due to no intensity in DNA FISH channel")
        else:
            threshold_min = otsu_threshold_val_local_img + (
                        img.max() - otsu_threshold_val_local_img) / min_threshold_factor

            img_seg_local = np.zeros_like(local_img_singlet)

            for k in range(local_cycle_factor):
                local = threshold_local(local_img, 10 * k + 7)
                out = (local_img_singlet > local).astype(int)
                out = binary_erosion(out)

                if filter == 'T':
                    if k == 0:
                        out = binary_dilation(out)
                        out_label = label(out)
                        out_props = regionprops(out_label, local_img_singlet)
                        for j in range(len(out_props)):
                            temp = np.zeros_like(local_img)
                            temp[out_label == out_props[j].label] = 1
                            temp_outer_edge = binary_dilation(temp, disk(6))
                            temp_outer_edge[temp == 1] = 0
                            mean_int_outer_edge = np.sum(local_img * temp_outer_edge) / np.sum(temp_outer_edge)
                            if (out_props[j].intensity_mean / mean_int_outer_edge > 1.2) & (out_props[j].area > min_size) & \
                                    (out_props[j].intensity_mean > min_threshold_first_round_factor * threshold_min):
                                img_seg_local[out_label == out_props[j].label] = 1
                    else:
                        out_label = label(out)
                        out_props = regionprops(out_label, local_img_singlet)
                        for j in range(len(out_props)):
                            if (out_props[j].intensity_mean > threshold_min) & (out_props[j].area > min_size):
                                img_seg_local[out_label == out_props[j].label] = 1
                elif filter == 'F':
                    img_seg_local[out == 1] = 1

            img_seg_watershed = np.zeros_like(local_img_singlet)
            bg_val = otsu_threshold_val_local_img * 3
            extreme_val = int(local_img_singlet.max() * 2 / otsu_threshold_val_local_img)
            maxima = extrema.h_maxima(local_img, extreme_val)
            elevation_map = sobel(local_img)
            markers = np.zeros_like(local_img)
            markers[local_img_singlet < bg_val] = 1
            markers[maxima == 1] = 2
            seg_wat = segmentation.watershed(elevation_map, markers)
            seg_wat_label = label(seg_wat)
            seg_wat_props = regionprops(seg_wat_label, local_img_singlet)
            for j in range(len(seg_wat_props)):
                if (seg_wat_props[j].intensity_mean > threshold_min) & (seg_wat_props[j].area > 12):
                    img_seg_watershed[seg_wat_label == seg_wat_props[j].label] = 1

            img_seg_individual = img_seg_watershed.copy()
            img_seg_individual[img_seg_local == 1] = 1
            img_seg_individual[local_nuclear_seg_convex == 0] = 0

            if training == 'T':
                # manuel correction
                viewer = napari.Viewer()
                viewer.add_image(local_img, blending='additive', colormap='green', contrast_limits=[0, img.max()])
                viewer.add_image(local_nuclear_seg_convex, blending='additive', colormap='blue')
                viewer.add_image(img_seg_individual, blending='additive', contrast_limits=[0, 1])
                shapes_signal_remove = viewer.add_shapes(name='signal to be removed', ndim=2)
                shapes_signal_add = viewer.add_shapes(name='signal to be added', ndim=2)
                napari.run()

                img_seg_individual = ima.napari_add_or_remove(shapes_signal_remove.data, 'remove', img_seg_individual)
                img_seg_individual = ima.napari_add_or_remove(shapes_signal_add.data, 'add', img_seg_individual)

            img_seg = ima.image_paste_to(img_seg, img_seg_individual,
                                         [int(original_centroid_nuclear[0] - local_centroid[0]),
                                          int(original_centroid_nuclear[1] - local_centroid[1])])
    return img_seg
