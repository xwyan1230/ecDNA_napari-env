import numpy as np
from skimage.morphology import remove_small_objects
from skimage.measure import regionprops, label
from matplotlib.path import Path


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


def points_exchange_xy(lst: list):
    """
    Exchange x and y for a list of points

    :param lst:
    :return:
    """
    out = []
    for i in range(len(lst)):
        out.append([lst[i][1], lst[i][0]])

    return out


def polygon_to_mask(img_shape: list, polygon: list):
    """
    Convert napari shape polygon into mask

    :param img_shape: image.shape()
    :param polygon: list of points
    :return:
    """
    poly_data = points_exchange_xy(polygon)
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


def blending(view):
    """
    For exporting napari view images (compress color images together)

    :param view: napari viewer
    :return:
    """
    blended = np.zeros(view.layers[0].data.shape + (4,))
    for layer in view.layers:
        # normalize data by clims
        normalized_data = (layer.data - layer.contrast_limits[0]) / (
                    layer.contrast_limits[1] - layer.contrast_limits[0])
        colormapped_data = layer.colormap.map(normalized_data.flatten())
        colormapped_data = colormapped_data.reshape(normalized_data.shape + (4,))

        blended = blended + colormapped_data

    blended[..., 3] = 1  # set alpha channel to 1
    blended[blended > 1] = 1
    blended[blended < 0] = 0
    return blended