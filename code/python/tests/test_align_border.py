#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# testAlignBorder.py: Test function for the border alignment alogrithm
# Author: Mathias Roesler
# Last modified: 11/23

import matplotlib.pyplot as plt
import numpy as np
import scipy.io

import thickness.projection as projection
import thickness.utils as utils

from thickness.constants import BASE, HOME


def alignBorderTest():
    """Tests the alignment alogrithm

    Args:

    Returns:

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

        for i, img in enumerate(img_stack):
            projection_points = projection.find_projection_points(
                img, centreline[i, :], params["nb_points"], params["horn"][i]
            )
            diff = np.diff(projection_points, axis=0)
            norm = np.linalg.norm(diff, axis=1)
            thickness = norm[np.arange(0, projection_points.shape[0], 2)]
            exclusion_indices = np.where(thickness == 0)[0]
            thickness[exclusion_indices] = float("nan")

            # Find the two halves
            indices = np.arange(len(projection_points))
            right_half = indices[(indices % 4 == 0) | ((indices - 1) % 4 == 0)]
            left_half = np.setdiff1d(indices, right_half)

            # Plot projection points first half
            plt.imshow(img, cmap="gray")
            plt.plot(
                projection_points[right_half][:, 0],
                projection_points[right_half][:, 1],
                ".b",
            )

            # Plot projection points second half
            plt.plot(
                projection_points[left_half][:, 0],
                projection_points[left_half][:, 1],
                ".g",
            )
            plt.show()

            # Order thickness to go from 0 to 2pi
            ordered_thickness = projection.align_border(thickness)

            x_values = np.linspace(
                0, 2 * np.pi, params["nb_points"], endpoint=False)
            plt.plot(x_values, ordered_thickness)
            plt.show()


alignBorderTest()
