#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# testLineCoordinates.py: Test function for the line coordinates algorithm
# Author: Mathias Roesler
# Last modified: 11/23

import matplotlib.pyplot as plt
import numpy as np
import scipy.io

import thickness_analysis.projection as projection
import utils.utils as utils


def lineCoordinatesTest():
    """Tests the line coordinate algorithm on test images

    Arguments:

    Return:

    """
    _dir = utils.HOME + "/" + utils.BASE + "/microCT/data/tests/"
    param_file = _dir + "test.toml"
    params = utils.parseTOML(param_file)
    cnt = 0  # Counter for theta
    _theta = [np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]

    for dataset in params["sets"]:
        print("Testing set {}".format(dataset))
        test_dir = _dir + dataset + "/muscle_segmentation"
        img_stack = utils.loadImageStack(test_dir)  # Load test images
        centreline_dict = scipy.io.loadmat(test_dir + "/centreline.mat")
        centreline = np.transpose(centreline_dict["centreline"])
        centreline = np.round(centreline).astype(int)  # Convert to int

        for i in range(2):
            if params["horn"][i] == "left":
                centre_point = centreline[i, 0:2]

            else:
                centre_point = centreline[i, 4:6]

            img = img_stack[i]
            line_x, line_y = projection.findLineCoordinates(
                img.shape, centre_point, _theta[cnt]
            )

            cnt += 1  # Increment counter

            plt.imshow(img, cmap="gray")
            plt.plot(line_x, line_y)
            plt.plot(centre_point[0], centre_point[1], ".r")
            plt.show()


lineCoordinatesTest()
