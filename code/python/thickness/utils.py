#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
utils.py

Utility functions for the thickness package
Author: Mathias Roesler
Date: 02/23
"""

import os
import sys
import glob
import tomli
import numpy as np
import skimage.io as skio


def load_image_stack(dir_path, extension="png"):
    """Loads the images found in the directory

    Args:
    dir_path -- str, path to the folder containing the images.
    extension -- str, image extension, default value png.

    Returns:
    img_stack -- ndarray, image stack.

    """
    if not os.path.isdir(dir_path):
        sys.stderr.write(
            "Error: the directory {} does not exists.\n".format(dir_path),
        )
        exit()

    img_list = sorted(glob.glob("*.{}".format(extension), root_dir=dir_path))

    # Load first image to get size of image
    img = skio.imread(os.path.join(dir_path, img_list[0]), as_gray=True)
    nb_x_pixels, nb_y_pixels = img.shape
    stack_size = len(img_list)

    # Pre-allocate and read images
    img_stack = np.zeros(
        [stack_size, nb_x_pixels, nb_y_pixels],
        dtype=np.uint8,
    )

    for i, img_name in enumerate(img_list):
        path = os.path.join(dir_path, img_name)

        if not os.path.isfile(path):
            sys.stderr.write("Error: {} is not a file.\n".format(path))
            exit()

        img_stack[i, :, :] = skio.imread(path, as_gray=True)

    return img_stack


def save_image_stack(
    img_stack,
    save_path,
    img_prefix,
    start_nb=0,
    extension="png",
):
    """Saves the images in the stack to the save directory

    Args:
    img_stack -- ndarray, stack of images to save.
    save_path -- str, path of the directory in which to save images.
    img_prefix -- str, prefix for the images in the stack.
    start_nb -- int, number at which to start saving the images,
            default value 0.
    extension -- str, image extension, default value png.

    Returns:

    """
    if not os.path.isdir(save_path):
        sys.stderr.write(
            "Error: the directory {} does not exists.\n".format(save_path),
        )
        exit()

    i = 0

    for img in img_stack:
        if not np.sum(img) == 0:
            img_path = os.path.join(
                save_path, img_prefix + str("%03d" % (i + start_nb))
            )

            if not img.dtype == np.dtype(np.uint8):
                img = img.astype(np.uint8)

            skio.imsave(
                "{}.{}".format(img_path, extension),
                img,
                check_contrast=False,
            )

            i += 1


def get_vector(p1, p2):
    """Finds the vector given two points

    Args:
    p1 -- np.array, first point.
    p2 -- np.array, second point.

    Returns:
    vec -- np.array, normalised vector between p1 and p2.

    """
    try:
        assert p1.shape == p2.shape

    except AssertionError:
        sys.stderr.write("Error: both points should have the same shape")
        exit()

    vec = p2.astype(np.int16) - p1.astype(np.int16)
    vec = vec / np.linalg.norm(vec)

    return vec


def get_angle(v1, v2):
    """Finds the angle given two vectors

    Args:
    v1 -- np.array, first vector.
    v2 -- np.array, second vector.

    Returns:
    angle -- float, angle between v1 and v2 in rad.

    """
    try:
        assert v1.shape == v2.shape

    except AssertionError:
        sys.stderr.write("Error: both vectors should have the same shape")
        exit()

    angle = np.arccos(
        np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)),
    )

    return angle


def parse_TOML(toml_file):
    """Parse a toml file

    Args:
    toml_file -- str, path to the toml file.

    Returns:
    data -- dict, dictonary containing the parameters.

    """
    with open(toml_file, "rb") as f:
        data = tomli.load(f)

    return data


def find_padding(cur_size, new_size):
    """Finds the padding size to add to be able to split an image

    Args:
    cur_size -- int, current size.
    new_size -- int, desired size.

    Returns:
    pad -- int, padding value to make the cur_size divisible by new_size.

    """
    if cur_size % new_size:
        return new_size * ((cur_size // new_size) + 1) - cur_size
    else:
        return 0


def moving_average(array, window_size):
    """Computes the moving average of an array

    Args:
    array -- ndarray, array over which to compute the average.
    window_size -- int, size of the window over which to compute.

    Returns:
    averaged_array -- ndarray, values of the moving average.

    """
    array_size = len(array)
    averaged_array = np.zeros(array.shape)
    half_window = window_size // 2

    if window_size > len(array):
        sys.stderr.write("Error: window size is greater than array size.\n")
        exit(1)

    elif window_size == len(array):
        return np.mean(array)

    for i in range(array_size):
        win_start = max(0, i - half_window)
        win_end = min(array_size, i + half_window + 1)

        window = array[win_start:win_end]
        averaged_array[i] = np.mean(window, axis=0)

    return averaged_array


def circular_average(array, window_size):
    """Computes the average for a circular array

    Args:
    array -- ndarray, array over which to compute the standard deviation.
    window_size -- int, size of the window over which to compute.

    Returns:
    mean_array -- ndarray, values of the moving standard deviation.

    """
    array_size = len(array)
    half_window_size = window_size // 2
    mean_array = np.zeros(array.shape)

    if window_size > len(array):
        sys.stderr.write("Error: window size is greater than array size.\n")
        exit(1)

    elif window_size == len(array):
        return [np.mean(array)]

    for i in range(len(array)):
        win_start = i - half_window_size
        win_end = i + half_window_size

        if abs(win_start) + abs(win_end) < window_size:
            # If window size is odd add 1
            win_end += 1

        # Calculate mean with wrapped edges
        mean_array[i] = np.mean(
            array[np.arange(win_start, win_end) % array_size], axis=0
        )

    return mean_array


def moving_std(array, window_size):
    """Computes the standard deviation for each window

    Args:
    array -- ndarray, array over which to compute the standard deviation.
    window_size -- int, size of the window over which to compute.

    Returns:
    std_array -- ndarray, values of the moving standard deviation.

    """
    array_size = len(array)
    std_array = np.zeros(array_size)

    if window_size > len(array):
        sys.stderr.write("Error: window size is greater than array size.\n")
        exit(1)

    elif window_size == len(array):
        return [np.std(array)]

    for i in range(0, array_size // window_size):
        window = array[i * window_size: (i + 1) * window_size]
        std_array[i * window_size + window_size // 5] = np.std(window)

    return std_array
