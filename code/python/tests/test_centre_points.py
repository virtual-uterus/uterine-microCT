#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
testCentrePoints.py

Test function for the centre point location
Author: Mathias Roesler
Date: 11/23
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.io

import thickness.utils as utils

from thickness.constants import BASE, HOME


def centrePointsTest():
    """Tests the location of the centre points in test images

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
            centre_point = centreline[i]
            img = img_stack[i]

            plt.imshow(img, cmap="gray")

            for j in range(3):
                plt.plot(centre_point[j * 2], centre_point[j * 2 + 1], ".r")

            plt.show()


centrePointsTest()
