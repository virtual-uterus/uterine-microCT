#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
testProjectionPoints.py

Test function for the projection points algorithm
Author: Mathias Roesler
Date: 11/23
"""

import numpy as np
import scipy.io

import thickness.plots as plots
import thickness.projection as projection
import thickness.utils as utils

from thickness.constants import BASE, HOME


def findProjectionPointsTest():
    """Tests the projection point algorithm on test images

    Arguments:

    Return:

    """
    _dir = HOME + "/" + BASE + "/uterine-microCT/data/tests/"
    param_file = _dir + "test.toml"
    params = utils.parse_TOML(param_file)

    for dataset in params["sets"]:
        print("Testing set {}".format(dataset))
        test_dir = _dir + dataset + "/muscle_segmentation"
        img_stack = utils.load_image_stack(test_dir)  # Load test images
        centreline_dict = scipy.io.loadmat(test_dir + "/centreline.mat")
        centreline = np.transpose(centreline_dict["centreline"])
        centreline = np.round(centreline).astype(int)  # Convert to int

        for i in range(2):
            projection_points = projection.find_projection_points(
                img_stack[i],
                centreline[i],
                params["nb_points"],
                params["horn"][i],
            )
            plots.plot_projection_points(
                img_stack[i], centreline[i, i *
                                         4: i * 4 + 2], projection_points
            )


findProjectionPointsTest()
