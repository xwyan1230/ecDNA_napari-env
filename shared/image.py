import numpy as np
from skimage.filters import threshold_otsu
from skimage.morphology import binary_dilation, binary_erosion, disk
import shared.objects as obj
import random
from matplotlib.path import Path
import imutils
import shared.dataframe as dat
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy
import napari


"""
# ---------------------------------------------------------------------------------------------------
# FUNCTIONS for IMAGE PROCESSING
# ---------------------------------------------------------------------------------------------------

img_to_int
    FUNCTION: convert RGB color image into intensity image
    SYNTAX:   img_to_int(img: np.array)

img_local_position
    FUNCTION: determine position of local image
    SYNTAX:   img_local_position(img: np.array, centroid: tuple, local_size: int)

img_local_seg
    FUNCTION: generate local segmented image
    SYNTAX:   img_local_seg(img: np.array, position: list, label_number: int)

edge_from_seg
    FUNCTION: detect edge from 0,1 image
    SYNTAX:   edge_from_seg(img: np.array)

image_paste
    FUNCTION: paste image based on reference distance to empty image
    SYNTAX:   image_paste(paste_to_size: np.array, paste_from_img: np.array, distance: list)

image_paste_to
    FUNCTION: paste image based on reference distance to a defined image
    SYNTAX:   image_paste_to(paste_to_image: np.array, paste_from_img: np.array, distance: list)

image_paste_fix_value
    FUNCTION: paste image based on reference distance for certain fixed value to a defined image
    SYNTAX:   image_paste_fix_value(paste_to_image: np.array, paste_from_img: np.array, distance: list, value: int)

image_deduction
    FUNCTION: image deduction (returns img1-img2, negative value replace as 0)
    SYNTAX:   image_deduction(img1: np.array, img2: np.array)

logical_ellipse
    FUNCTION: draw logical ellipse on given image
    SYNTAX:   logical_ellipse(img: np.array, centerX: int, centerY: int, a: int, b: int, avg=1)

distance_map_from_point
    FUNCTION: generate distance map from point for given image
    SYNTAX:   distance_map_from_point(img: np.array, point: tuple)

radial_distribution_from_distance_map
    FUNCTION: calculate radial distribution from distance map
    SYNTAX:   radial_distribution_from_distance_map(img_seg: np.array, img_distance_map: np.array, img_feature: 
              np.array, interval: float, feature_max: float)

save_3d_image
    FUNCTION: save 3d images as txt file (require reshape)
    SYNTAX:   save_3d_image(img_stack: np.array, save_location: str, file_name: str)

load_3d_image
    FUNCTION: load 3d image as np.array (require reshape)
    SYNTAX:   load_3d_image(load_location: str, file_name: str, shape2: int)

angle_map_from_point
    FUNCTION: generate angle map from point for given image
    SYNTAX:   angle_map_from_point(img: np.array, point: tuple)

polygon_to_mask
    FUNCTION: convert napari shape polygon into mask
    SYNTAX:   polygon_to_mask(img_shape: list, polygon: list)

napari_add_or_remove
    FUNCTION: manual correction by adding or removing objects based on napari shapes
    SYNTAX:   napari_add_or_remove(shape_data, option: str, modify_mask: np.array)
    
napari_add_or_remove_obj
    FUNCTION: manual correction by adding or removing objects based on napari shapes
    SYNTAX:   manual correction by adding or removing objects based on napari shapes

napari_change_between_masks
    FUNCTION: manual correction by moving objects between masks based on napari shapes
    SYNTAX:   napari_change_between_masks(shape_data, change_from_mask: np.array, change_to_mask: np.array)

img_background_correction
    FUNCTION: generate background corrected image based on background intensity subtraction
    SYNTAX:   img_background_correction(img: np.array, bg_int: float)
"""


def img_to_int(img: np.array):
    """
    Convert RGB color image into intensity image

    :param img: color image acquired from EVOS-M5000, could be any channel including TRANS
    :return: img: intensity image
    """
    out = img
    if len(np.shape(img)) == 3:
        out = img[:, :, 0] + img[:, :, 1] + img[:, :, 2]
    elif len(np.shape(img)) == 4:
        out = img[:, :, :, 0] + img[:, :, :, 1] + img[:, :, :, 2]
    return out


def img_local_position(img: np.array, centroid: tuple, local_size: int):
    """
    Determine position of local image

    Check if local image is within the original image, otherwise, substitute the boundary with original image

    :param img: np.array, original image
    :param centroid: tuple, center of the new local image
    :param local_size: int, half size of the new local image
    :return: return the four boundary of the new image
    """
    shape = img.shape
    left = int(centroid[0] - local_size) if centroid[0] > local_size else 0
    right = int(centroid[0] + local_size) if centroid[0] + local_size < shape[0] else shape[0]
    top = int(centroid[1] - local_size) if centroid[1] > local_size else 0
    bottom = int(centroid[1] + local_size) if centroid[1] + local_size < shape[1] else shape[1]

    return left, right, top, bottom


def img_local_seg(img: np.array, position: list, label_number: int):
    """
    Generate local segmented image

    :param img: np.array, original segmented image
    :param position: list, position of local segmented image
    :param label_number: int, number of label in original segmented image
    :return:
    """
    temp = img.copy()
    temp = temp[position[0]:position[1], position[2]:position[3]]
    out = np.zeros_like(temp, dtype=float)
    out[temp == label_number] = 1

    return out


def img_local_seg_center(img: np.array, position: list):
    """
    Generate local segmented image

    :param img: np.array, original segmented image
    :param position: list, position of local segmented image
    :return:
    """
    temp = img.copy()
    temp = temp[position[0]:position[1], position[2]:position[3]]
    out = np.zeros_like(temp, dtype=float)
    label_number = temp[int(out.shape[0]/2)][int(out.shape[1]/2)]
    out[temp == label_number] = 1

    return out


def edge_from_seg(img: np.array):
    """
    Detect edge from 0,1 image

    Purpose: the boundary of region 1 will be labeled as 1

    :param img: np.array, original 0, 1 image
    :return:
    """
    shape = img.shape
    out = np.zeros_like(img, dtype=float)
    for m in range(shape[0]-2):
        for n in range(shape[1]-2):
            if (img[m+1, n+1] == 1) & ((img[m, n] == 0)|(img[m, n+1] == 0)|(img[m, n+2] == 0)|(img[m+1, n] == 0)|
                                       (img[m+1, n+2] == 0)|(img[m+2, n] == 0)|(img[m+2, n+1] == 0)|(img[m+2, n+2] == 0)):
                out[m+1, n+1] = 1

    return out


def image_paste(paste_to_size: np.array, paste_from_img: np.array, distance: list):
    """
    Paste image based on reference distance to empty image

    :param paste_to_size: np.array, same size image as output
    :param paste_from_img: np.array, image to be pasted
    :param distance: [x,y] reference distance
    :return:
    """
    paste_to_img = np.zeros_like(paste_to_size)
    for i in range(paste_from_img.shape[0]):
        for j in range(paste_from_img.shape[1]):
            if paste_from_img[i][j] != 0:
                paste_to_img[i+distance[0]][j+distance[1]] = paste_from_img[i][j]

    return paste_to_img


def image_paste_to(paste_to_image: np.array, paste_from_img: np.array, distance: list):
    """
    Paste image based on reference distance to a defined image

    :param paste_to_image: np.array, image gets pasted to
    :param paste_from_img: np.array, image to be pasted
    :param distance: [x,y] reference distance
    :return:
    """
    paste_to_img = paste_to_image.copy()
    for i in range(paste_from_img.shape[0]):
        for j in range(paste_from_img.shape[1]):
            if paste_from_img[i][j] != 0:
                paste_to_img[i+distance[0]][j+distance[1]] = paste_from_img[i][j]

    return paste_to_img


def image_paste_fix_value(paste_to_image: np.array, paste_from_img: np.array, distance: list, value: int):
    """
        Paste image based on reference distance for certain fixed value to a defined image

        :param paste_to_image: np.array, image gets pasted to
        :param paste_from_img: np.array, image to be pasted (only shape is pasted)
        :param distance: [x,y] reference distance
        :param value: fixed value
        :return:
        """
    paste_to_img = paste_to_image.copy()
    for i in range(paste_from_img.shape[0]):
        for j in range(paste_from_img.shape[1]):
            if (paste_from_img[i][j] != 0) & (i+distance[0] < len(paste_to_img)) & \
                    (j+distance[1] < len(paste_to_img[0])):
                paste_to_img[i + distance[0]][j + distance[1]] = value

    return paste_to_img


def image_deduction(img1: np.array, img2: np.array):
    """
    Image deduction (returns img1-img2, negative value replace as 0)

    :param img1: np.array
    :param img2: np.array
    :return:
    """
    out = np.zeros_like(img1)
    for i in range(len(img1)):
        for j in range(len(img1[i])):
            if img1[i][j] > img2[i][j]:
                out[i][j] = img1[i][j] - img2[i][j]

    return out


def logical_ellipse(img: np.array, centerX: int, centerY: int, a: int, b: int, avg=1):
    """
    Draw logical ellipse on given image

    :param img: np.array, input image
    :param centerX: center coordinate X
    :param centerY: center coordinate Y
    :param a: r on X axis
    :param b: r on Y axis
    :param avg: intensity value for the ellipse, default=1
    :return:
    """
    out = img.copy()
    for i in range(len(img)):
        for j in range(len(img[0])):
            if (i-centerX)**2/a**2 + (j-centerY)**2/b**2 <= 1:
                out[i][j] = avg

    return out


def logical_dot_sample(img: np.array, mask: np.array, n: int, avg=1):
    """
    Create n randomly selected dot within certain region

    :param img: np.array, input image
    :param mask: np.array, mask image to determine the region
    :param n: int, number of dots
    :param avg: intensity value for the ellipse, default=1
    :return:
    """
    out = img.copy()
    mask_int = np.sum(mask)
    temp_lst = list(np.arange(mask_int))
    target_lst = random.sample(temp_lst, n)
    k = 0
    for i in range(len(img)):
        for j in range(len(img[0])):
            if mask[i][j] == 1:
                if k in target_lst:
                    out[i][j] = avg
                k = k+1

    return out


def distance_map_from_point(img: np.array, point: tuple):
    """
    Generate distance map from point for given image

    :param img: np.array, input image
    :param point: given point
    :return:
    """
    p = [round(point[0]), round(point[1])]
    image_shape = [len(img), len(img[0])]
    ydis = np.array([np.arange(-p[0], image_shape[0] - p[0], 1)]).transpose() * np.ones(image_shape[1])
    xdis = np.array([np.ones(image_shape[0])]).transpose() * np.arange(-p[1], image_shape[1] - p[1], 1)
    out = (ydis * ydis + xdis * xdis)**0.5

    return out


def radial_distribution_from_distance_map(img_seg: np.array, img_distance_map: np.array, img_feature: np.array,
                                          interval: float, feature_max: float):
    """
    Calculate radial distribution from distance map

    :param img_seg: np.array, 0-1 binary image, segmentation image
    :param img_distance_map: np.array, distance map
    :param img_feature: np.array, feature image, for example: intensity image
    :param interval: float, bin size
    :param feature_max: float, maximum for analysis, number larger than maximum number will be binned in the last bin
    :return:
    """
    out = []
    """feature_seg = np.array(img_feature.copy()).astype(int)
    feature_seg = feature_seg - feature_thresh
    feature_seg[(feature_seg < 0) | (img_seg == 0)] = 0
    feature_thresh = img_feature.copy()
    feature_thresh[feature_seg == 0] = 0
    seg_seg = img_seg.copy()
    seg_seg[feature_seg == 0] = 0"""
    feature_thresh = img_feature.copy()
    feature_thresh[img_seg == 0] = 0
    seg_seg = img_seg.copy()
    mean_feature = np.sum(feature_thresh)/np.sum(seg_seg)
    for i in np.arange(0, feature_max, interval):
        seg_temp = seg_seg.copy()
        feature_temp = feature_thresh.copy()
        seg_temp[img_distance_map < i] = 0
        feature_temp[img_distance_map < i] = 0
        if feature_max - i > interval:
            seg_temp[img_distance_map >= (i+interval)] = 0
            feature_temp[img_distance_map >= (i+interval)] = 0
        if np.sum(seg_temp) != 0:
            out.append((np.sum(feature_temp)/np.sum(seg_temp))/mean_feature)
        else:
            out.append(0)

    return out


def sum_up_image(img_add_to: np.array, img_add_from: np.array, direction: tuple, normalization_factor: float):
    """
    Sum up image (need to rescale before display)

    :param img_add_to: np.array
    :param img_add_from: np.array
    :param direction: distance between two anchor points of two images
    :param normalization_factor: float
    :return:
    """
    feature_recentered = image_paste(img_add_to, img_add_from, direction)
    feature_recentered = np.array(feature_recentered).astype(float) * normalization_factor
    out = img_add_to + feature_recentered

    return out


def radial_percentage_from_distance_map(img_seg: np.array, img_distance_map: np.array, img_feature: np.array,
                                        feature_range: list):
    """
    Get radial percentage within certain range from distance map

    :param img_seg: np.array, 0-1 binary image, segmentation image
    :param img_distance_map: np.array, distance map
    :param img_feature: np.array, feature image, for example: intensity image
    :param feature_range: list of two values, [lower limit, higher limit)
    :return:
    """
    feature_thresh = img_feature.copy()
    feature_thresh[img_seg == 0] = 0
    sum_feature = np.sum(feature_thresh)

    feature_temp = feature_thresh.copy()
    feature_temp[(img_distance_map < feature_range[0]) | (img_distance_map >= feature_range[1])] = 0

    return np.sum(feature_temp)*1.0/sum_feature


def save_3d_image(img_stack: np.array, save_location: str, file_name: str):
    """
    Save 3d images as txt file (require reshape)
    :param img_stack: np.array, 3d image stack, time or z
    :param save_location: str, master_folder
    :param file_name: str, saving name
    :return:
    """
    arr_reshaped = img_stack.reshape(img_stack.shape[0], -1)
    np.savetxt("%s%s.txt" % (save_location, file_name), arr_reshaped)


def load_3d_image(load_location: str, file_name: str, shape2: int):
    """
    Load 3d image as np.array (require reshape)
    :param load_location: str, master_folder
    :param file_name: str, saving name
    :param shape2: int, original_array.shape[2]
    :return:
    """
    loaded_arr = np.loadtxt("%s%s.txt" % (load_location, file_name))
    load_original_arr = loaded_arr.reshape(loaded_arr.shape[0], loaded_arr.shape[1] // shape2, shape2)
    return load_original_arr


def angle_map_from_point(img: np.array, point: tuple):
    """
    Generate angle map from point for given image

    :param img: np.array, input image
    :param point: given point
    :return:
    """
    p = [round(point[0]), round(point[1])]
    image_shape = [len(img), len(img[0])]
    ydis = np.array([np.arange(-p[0], image_shape[0] - p[0], 1)]).transpose() * np.ones(image_shape[1])
    xdis = np.array([np.ones(image_shape[0])]).transpose() * np.arange(-p[1], image_shape[1] - p[1], 1)
    dis = (ydis * ydis + xdis * xdis)**0.5
    np.seterr(invalid='ignore')
    out = np.arcsin(np.abs(xdis)/dis) * 180/np.pi
    out[(ydis > 0) & (xdis >= 0)] = 180 - out[(ydis > 0) & (xdis >= 0)]
    out[(ydis > 0) & (xdis < 0)] = out[(ydis > 0) & (xdis < 0)] + 180
    out[(ydis <= 0) & (xdis < 0)] = 360 - out[(ydis <= 0) & (xdis < 0)]

    return out


def polygon_to_mask(img_shape: list, polygon: list):
    """
    Convert napari shape polygon into mask

    :param img_shape: image.shape()
    :param polygon: list of points
    :return:
    """
    poly_data = dat.points_exchange_xy(polygon)
    nx = img_shape[1]
    ny = img_shape[0]
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))
    x, y = x.flatten(), y.flatten()
    points = np.vstack((x, y)).T
    path = Path(poly_data)
    grid = path.contains_points(points)
    grid = grid.reshape(ny, nx)
    mask = grid

    return mask


def napari_add_or_remove(shape_data, option: str, modify_mask: np.array):
    """
    Manual correction by adding or removing objects based on napari shapes

    :param shape_data: polygon data from napari shape
    :param option: only accept 'add' or 'remove'
    :param modify_mask: the mask that is to be modified
    :return:
    """
    out = modify_mask.copy()
    for i in range(len(shape_data)):
        poly_data = shape_data[i]
        mask = polygon_to_mask(modify_mask.shape, poly_data)
        if option == 'add':
            out[mask == 1] = 1
        elif option == 'remove':
            out[mask == 1] = 0

    return out


def napari_add_or_remove_obj(shape_data, option: str, modify_mask: np.array):
    """
    Manual correction by adding or removing objects based on napari shapes

    :param shape_data: polygon data from napari shape
    :param option: only accept 'add' or 'remove'
    :param modify_mask: the mask that is to be modified, labeled image
    :return:
    """
    out = modify_mask.copy()
    max_out = modify_mask.max()
    for i in range(len(shape_data)):
        poly_data = shape_data[i]
        mask = polygon_to_mask(modify_mask.shape, poly_data)
        if option == 'add':
            out[mask == 1] = max_out+1+i
        elif option == 'remove':
            out[mask == 1] = 0

    return out


def napari_change_between_masks(shape_data, change_from_mask: np.array, change_to_mask: np.array):
    """
    Manual correction by moving objects between masks based on napari shapes

    :param shape_data: polygon data from napari shape
    :param change_from_mask: the mask that is to be removed
    :param change_to_mask: the mask that is to be added
    :return:
    """
    out_from = change_from_mask.copy()
    out_to = change_to_mask.copy()
    for i in range(len(shape_data)):
        poly_data = shape_data[i]
        mask = polygon_to_mask(change_from_mask.shape, poly_data)
        out_to[mask * change_from_mask == 1] = 1
        out_from[mask == 1] = 0

    return out_from, out_to


def img_background_correction(img: np.array, bg_int: float):
    """
    Generate background corrected image based on background intensity subtraction

    :param img: np.array
    :param bg_int: float
    :return:
    """
    out = img.copy().astype(float) - np.ones_like(img) * bg_int
    out[out < 0] = 0
    out = out.astype(int)

    return out


def img_align_move(img1: np.array, img2: np.array, rotate_a, topleft_target):
    img1_out = img1.copy()
    img1_temp = imutils.rotate(img1_out, angle=rotate_a[0])
    img1_out = img1_temp
    img2_out = img2.copy()
    img2_temp = imutils.rotate(img2_out, angle=rotate_a[1])
    img2_out = img2_temp
    if topleft_target[0] > 0:
        img2_temp = np.concatenate([np.zeros(shape=[img2_out.shape[0], abs(int(topleft_target[0]))]), img2_out], axis=1)
        img2_out = img2_temp
    else:
        img1_temp = np.concatenate([np.zeros(shape=[img1_out.shape[0], abs(int(topleft_target[0]))]), img1_out], axis=1)
        img1_out = img1_temp
    if topleft_target[1] > 0:
        img2_temp = np.concatenate([np.zeros(shape=[abs(int(topleft_target[1])), img2_out.shape[1]]), img2_out], axis=0)
        img2_out = img2_temp
    else:
        img1_temp = np.concatenate([np.zeros(shape=[abs(int(topleft_target[1])), img1_out.shape[1]]), img1_out], axis=0)
        img1_out = img1_temp
    img1_out, img2_out = img_same_size(img1_out, img2_out)
    return img1_out, img2_out


def img_search_global(img: np.array, img_search: np.array, poly_data, interval):
    topleft_ori = [poly_data[0][1], poly_data[0][0]]
    xrange = poly_data[1][1] - poly_data[0][1]
    yrange = poly_data[2][0] - poly_data[0][0]

    min_ratio = 1
    topleft_target = topleft_ori
    i = 0
    for x in range(int(xrange / interval)):
        for y in range(int(yrange / interval)):
            i = i + 1
            topleft = [topleft_ori[0] + x * interval, topleft_ori[1] + y * interval]
            img_cut = img.copy()
            img_cut = img_cut[int(topleft[1]):int(topleft[1] + img_search.shape[0]),
                      int(topleft[0]):int(topleft[0] + img_search.shape[1])]

            minus = img_search.astype(float) - img_cut.astype(float)
            minus[minus < 0] = 0
            min_ratio_temp = minus.sum() / img_search.sum()
            print('%s/%s: %s' % (i, int(xrange / interval) * int(yrange / interval), min_ratio_temp))
            if min_ratio_temp < min_ratio:
                min_ratio = min_ratio_temp
                topleft_target = [topleft_ori[0] + x * interval, topleft_ori[1] + y * interval]
    print(min_ratio)
    print(topleft_target)
    return topleft_target, min_ratio


def img_search_local(img: np.array, img_search: np.array, topleft_center, interval):
    topleft_ori = [topleft_center[0] - 10*interval, topleft_center[1]-10*interval]
    xrange = 20 * interval
    yrange = 20 * interval

    min_ratio = 1
    topleft_target = topleft_ori
    i = 0
    for x in range(int(xrange / interval)):
        for y in range(int(yrange / interval)):
            i = i + 1
            topleft = [topleft_ori[0] + x * interval, topleft_ori[1] + y * interval]
            img_cut = img.copy()
            img_cut = img_cut[int(topleft[1]):int(topleft[1] + img_search.shape[0]),
                      int(topleft[0]):int(topleft[0] + img_search.shape[1])]

            minus = img_search.astype(float) - img_cut.astype(float)
            minus[minus < 0] = 0
            min_ratio_temp = minus.sum() / img_search.sum()
            print('%s/%s: %s' % (i, int(xrange / interval) * int(yrange / interval), min_ratio_temp))
            if min_ratio_temp < min_ratio:
                min_ratio = min_ratio_temp
                topleft_target = [topleft_ori[0] + x * interval, topleft_ori[1] + y * interval]
    print(min_ratio)
    print(topleft_target)
    return topleft_target, min_ratio


def img_same_size(img1: np.array, img2: np.array):
    shape_diff = [img1.shape[0] - img2.shape[0], img1.shape[1] - img2.shape[1]]
    img2_out = img2.copy()
    img1_out = img1.copy()
    if shape_diff[0] > 0:
        img2_temp = np.concatenate([img2_out, np.zeros(shape=[abs(int(shape_diff[0])), img2_out.shape[1]])], axis=0)
        img2_out = img2_temp
    else:
        img1_temp = np.concatenate([img1_out, np.zeros(shape=[abs(int(shape_diff[0])), img1_out.shape[1]])], axis=0)
        img1_out = img1_temp
    if shape_diff[1] > 0:
        img2_temp = np.concatenate([img2_out, np.zeros(shape=[img2_out.shape[0], abs(int(shape_diff[1]))])], axis=1)
        img2_out = img2_temp
    else:
        img1_temp = np.concatenate([img1_out, np.zeros(shape=[img1_out.shape[0], abs(int(shape_diff[1]))])], axis=1)
        img1_out = img1_temp
    return img1_out, img2_out


def img_search_global10(img: np.array, img_search: np.array, poly_data, interval):
    topleft_ori = [poly_data[0][1], poly_data[0][0]]
    xrange = poly_data[1][1] - poly_data[0][1]
    yrange = poly_data[2][0] - poly_data[0][0]

    min_ratio = []
    topleft_target = []
    i = 0
    for x in range(int(xrange / interval)):
        for y in range(int(yrange / interval)):
            i = i + 1
            topleft = [topleft_ori[0] + x * interval, topleft_ori[1] + y * interval]
            img_cut = img.copy()
            img_cut = img_cut[int(topleft[1]):int(topleft[1] + img_search.shape[0]),
                      int(topleft[0]):int(topleft[0] + img_search.shape[1])]

            minus = img_search.astype(float) - img_cut.astype(float)
            minus[minus < 0] = 0
            min_ratio_temp = minus.sum() / img_search.sum()
            print('%s/%s: %s' % (i, int(xrange / interval) * int(yrange / interval), min_ratio_temp))
            min_ratio.append(min_ratio_temp)
            topleft_target.append([topleft_ori[0] + x * interval, topleft_ori[1] + y * interval])
    return topleft_target[:10], min_ratio[:10]