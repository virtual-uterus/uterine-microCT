#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# utils.py: Utility functions for the quantative analysis
# Author: Mathias Roesler
# Last modified: 02/23

import os
import sys
import glob
import tomli
import numpy as np
import skimage.io as skio

HOME = os.path.expanduser('~')
BASE = "Documents/phd"


def loadImageStack(dir_path, extension="png"):
    """ Loads the images found in the directory

    Arguments:
    dir_path -- str, path to the folder containing the images.
    extension -- str, image extension, default value png.

    Return:
    img_stack -- ndarray, image stack.

    """
    if not os.path.isdir(dir_path):
        sys.stderr.write("Error: the directory {} does not exists.\n".format(
            dir_path))
        exit()

    img_list = sorted(glob.glob("*.{}".format(extension), root_dir=dir_path))

    # Load first image to get size of image
    img = skio.imread(os.path.join(dir_path, img_list[0]), as_gray=True)
    nb_x_pixels, nb_y_pixels = img.shape
    stack_size = len(img_list)

    # Pre-allocate and read images
    img_stack = np.zeros(
        [stack_size, nb_x_pixels, nb_y_pixels], dtype=np.uint8)

    for i, img_name in enumerate(img_list):
        path = os.path.join(dir_path, img_name)

        if not os.path.isfile(path):
            sys.stderr.write("Error: {} is not a file.\n".format(path))
            exit()

        img_stack[i, :, :] = skio.imread(path, as_gray=True)

    return img_stack


def saveImageStack(img_stack, save_path, img_prefix, start_nb=0,
                   extension="png"):
    """ Saves the images in the stack to the save directory

    Arguments:
    img_stack -- ndarray, stack of images to save.
    save_path -- str, path of the directory in which to save images.
    img_prefix -- str, prefix for the images in the stack.
    start_nb -- int, number at which to start saving the images,
            default value 0.
    extension -- str, image extension, default value png.

    Return:

    """
    if not os.path.isdir(save_path):
        sys.stderr.write("Error: the directory {} does not exists.\n".format(
            save_path))
        exit()

    i = 0

    for img in img_stack:
        if not np.sum(img) == 0:
            img_path = os.path.join(save_path, img_prefix + str(
                "%03d" % (i + start_nb)))

            if not img.dtype == np.dtype(np.uint8):
                img = img.astype(np.uint8)

            skio.imsave("{}.{}".format(img_path, extension), img,
                        check_contrast=False)

            i += 1


def getVector(p1, p2):
    """ Finds the vector given two points

    Arguments:
    p1 -- np.array, first point.
    p2 -- np.array, second point.

    Return:
    vec -- np.array, normalised vector between p1 and p2.

    """
    try:
        assert (p1.shape == p2.shape)

    except AssertionError:
        sys.stderr.write("Error: both points should have the same shape")
        exit()

    vec = p2.astype(np.int16) - p1.astype(np.int16)
    vec = vec / np.linalg.norm(vec)

    return vec


def getAngle(v1, v2):
    """ Finds the angle given two vectors

    Arguments:
    v1 -- np.array, first vector.
    v2 -- np.array, second vector.

    Return:
    angle -- float, angle between v1 and v2 in rad.

    """
    try:
        assert (v1.shape == v2.shape)

    except AssertionError:
        sys.stderr.write("Error: both vectors should have the same shape")
        exit()

    angle = np.arccos(np.dot(v1, v2) / (
        np.linalg.norm(v1) * np.linalg.norm(v2))
    )

    return angle


def parseTOML(toml_file):
    """ Parse a toml file

    Arguments:
    toml_file -- str, path to the toml file.

    Return:
    data -- dict, dictonary containing the parameters.

    """
    with open(toml_file, 'rb') as f:
        data = tomli.load(f)

    return data


def findPadding(cur_size, new_size):
    """ Finds the padding size to add to be able to split an image

    Arguments:
    cur_size -- int, current size.
    new_size -- int, desired size.

    Return:
    pad -- int, padding value to make the cur_size divisible by new_size.

    """
    if cur_size % new_size:
        return new_size * ((cur_size // new_size) + 1) - cur_size
    else:
        return 0


def movingAverage(array, window_size):
    """ Computes the moving average of an array

    Arguments:
    array -- ndarray, array over which to compute the average.
    window_size -- int, size of the window over which to compute.

    Return:
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
        win_start = max(0, i-half_window)
        win_end = min(array_size, i+half_window+1)

        window = array[win_start:win_end]
        averaged_array[i] = np.mean(window, axis=0)

    return averaged_array


def circularAverage(array, window_size):
    """ Computes the average for a circular array

    Arguments:
    array -- ndarray, array over which to compute the standard deviation.
    window_size -- int, size of the window over which to compute.

    Return:
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
        mean_array[i] = np.mean(array[
            np.arange(win_start, win_end) % array_size], axis=0)

    return mean_array


def movingStd(array, window_size):
    """ Computes the standard deviation for each window

    Arguments:
    array -- ndarray, array over which to compute the standard deviation.
    window_size -- int, size of the window over which to compute.

    Return:
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
        window = array[i*window_size:(i+1)*window_size]
        std_array[i*window_size + window_size // 5] = np.std(window)

    return std_array


def writeExElemVol(file_path, elements, thickness=True):
    """ Writes out the data from a volumetric mesh to a exnode file

    Arguments:
    file_path -- str, path to the file to save to.
    elements -- ndarray, list of nodes associated with each tetrahedra,
            size = Nx4.
    thickness -- bool, flag used if thickness has been provided to the
            exnode file, default True.

    Return:

    """
    try:
        assert (elements.shape[1] == 4)

    except AssertionError:
        sys.stderr.write("Error: elements should contain 4 nodes\n")

    with open(file_path, "w") as f:
        # Write the exnode file header
        f.write("Group name: mesh\n")
        f.write("Region: /uterus\n")
        f.write("Shape.  Dimension=3 simplex(2;3)*simplex*simplex\n")
        f.write("#Scale factor sets=0\n")
        f.write("#Nodes=4\n")

        if thickness:
            f.write("#Fields=2\n")

        else:
            # If no thickness is provided there is only one field
            f.write("#Fields=1\n")

        f.write(
            "1) coordinates, coordinate, rectangular cartesian, "
            "#Components=3\n")
        f.write(
            " x. l.simplex(2;3)*l.simplex*l.simplex, no modify, "
            "standard node based.\n")
        f.write("  #Nodes=4\n")
        f.write("  1. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write("  2. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write("  3. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write("  4. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write(
            " y. l.simplex(2;3)*l.simplex*l.simplex, no modify, "
            "standard node based.\n")
        f.write("  #Nodes=4\n")
        f.write("  1. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write("  2. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write("  3. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write("  4. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write(
            " z. l.simplex(2;3)*l.simplex*l.simplex, no modify, "
            "standard node based.\n")
        f.write("  #Nodes=4\n")
        f.write("  1. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write("  2. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write("  3. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write("  4. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")

        if thickness:
            f.write(
                "2) thickness, field, rectangular cartesian, "
                "#Components=1\n")
            f.write(" thickness. constant, no modify, standard node based.\n")
            f.write("  #Nodes=1\n")
            f.write("  1. #Values=1\n")
            f.write("	Value indices: 1\n")
            f.write("	Scale factor indices: 0\n")

        for i, nodes in enumerate(elements):
            f.write("Element: {} 0 0\n".format(i+1))
            f.write(" Nodes: \n")
            f.write("  {} {} {} {}\n".format(
                    nodes[0]+1, nodes[1]+1, nodes[2]+1, nodes[3]+1))


def writeExElemSurf(file_path, elements, thickness=True):
    """ Writes out the data from a surface mesh to a exnode file

    Arguments:
    file_path -- str, path to the file to save to.
    elements -- ndarray, list of nodes associated with each triangle,
            size = Nx3.
    thickness -- bool, flag used if thickness has been provided to the
            exnode file, default True.

    Return:

    """
    try:
        assert (elements.shape[1] == 3)

    except AssertionError:
        sys.stderr.write("Error: elements should contain 3 nodes\n")

    with open(file_path, "w") as f:
        # Write the exnode file header
        f.write("Group name: mesh\n")
        f.write("Region: /uterus\n")
        f.write("Shape.  Dimension=2 simplex(2)*simplex\n")
        f.write("#Scale factor sets=0\n")
        f.write("#Nodes=3\n")

        if thickness:
            f.write("#Fields=2\n")

        else:
            # If no thickness is provided there is only one field
            f.write("#Fields=1\n")

        f.write(
            "1) coordinates, coordinate, rectangular cartesian, "
            "#Components=3\n")
        f.write(
            " x. l.simplex(2)*l.simplex, no modify, standard node based.\n")
        f.write("  #Nodes=3\n")
        f.write("  1. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write("  2. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write("  3. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write(
            " y. l.simplex(2)*l.simplex, no modify, standard node based.\n")
        f.write("  #Nodes=3\n")
        f.write("  1. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write("  2. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write("  3. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write(
            " z. l.simplex(2)*l.simplex, no modify, standard node based.\n")
        f.write("  #Nodes=3\n")
        f.write("  1. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write("  2. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")
        f.write("  3. #Values=1\n")
        f.write("	Value indices: 1\n")
        f.write("	Scale factor indices: 0\n")

        if thickness:
            f.write(
                "2) thickness, field, rectangular cartesian, #Components=1\n")
            f.write(" thickness. constant, no modify, standard node based.\n")
            f.write("  #Nodes=1\n")
            f.write("  1. #Values=1\n")
            f.write("	Value indices: 1\n")
            f.write("	Scale factor indices: 0\n")

        for i, nodes in enumerate(elements):
            f.write("Element: {} 0 0\n".format(i+1))
            f.write(" Nodes: \n")
            f.write("  {} {} {}\n".format(
                    nodes[0]+1, nodes[1]+1, nodes[2]+1))


def writeExNode(file_path, nodes, thickness=None):
    """ Writes out the nodes from a mesh to a exnode file,
            and adds the thickness field if provided

    Arguments:
    file_path -- str, path to the file to save to.
    nodes -- ndarray, list of coordinates for each node.
            size = Nx3
    thickness -- ndarray, list of thickness value for each node.
            size = Nx1, default value None.

    Return:

    """
    try:
        # Check for number of coordinates
        assert (nodes.shape[1] == 3)

    except AssertionError:
        sys.stderr.write("Error: nodes should have three coordinates\n")
        exit()

    if type(thickness) is not type(None):
        try:
            # Check that thickness and nodes have the same dimension
            assert (nodes.shape[0] == thickness.shape[0])

        except AssertionError:
            sys.stderr.write("Error: nodes and thickness should have the same "
                             "number of elements\n")
            exit()

    with open(file_path, "w") as f:
        # Write exnode file header
        f.write("Group name: mesh\n")
        f.write("Region: /uterus\n")

        if type(thickness) is not type(None):
            f.write("#Fields=2\n")

        else:
            # If no thickness is provided there is only one field
            f.write("#Fields=1\n")

        f.write(
            "1) coordinates, coordinate, rectangular cartesian, "
            "#Components=3\n")
        f.write(" x. Value index=1, #Derivatives=0\n")
        f.write(" y. Value index=2, #Derivatives=0\n")
        f.write(" z. Value index=3, #Derivatives=0\n")

        if type(thickness) is not type(None):
            f.write(
                "2) thickness, field, rectangular cartesian, #Components=1\n")
            f.write(" thickness. Value index=4, #Derivatives=0\n")

        for i in range(len(nodes)):
            f.write("Node: {}\n".format(i+1))
            f.write(" {} {} {}\n".format(
                    nodes[i][0], nodes[i][1], nodes[i][2]))

            if type(thickness) is not type(None):
                f.write(" {}\n".format(thickness[i]))
