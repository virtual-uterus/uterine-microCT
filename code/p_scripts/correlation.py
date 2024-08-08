#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# correlation.py: Script to estimate correlation between histology and uCT
# Author: Mathias Roesler
# Last modified: 06/23

import argparse
import os
import sys

import numpy as np
from scipy import stats

import utils.utils as utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Specific script to estimate the correlation between "
        " histology and uCT"
    )

    parser.add_argument(
        "uCT_path", type=str, metavar="uCT-path",
        help="path from BASE to the uCT data"
    )
    parser.add_argument(
        "histo_path",
        type=str,
        metavar="histo-path",
        help="path from BASE to the histology data",
    )
    parser.add_argument(
        "base_name", type=str, metavar="base-name", help="name of the dataset"
    )
    parser.add_argument(
        "--horn",
        type=str,
        choices={"left", "right"},
        help="horn to process, default right",
        default="right",
    )
    parser.add_argument(
        "--not-d",
        action="store_true",
        help="flag used if the uCT dataset is not downsampled, default False",
    )

    # Parse input arguments
    args = parser.parse_args()

    # Set up variables
    uCT_directory = os.path.join(utils.HOME, utils.BASE, args.uCT_path,
                                 args.base_name)
    histo_directory = os.path.join(
        utils.HOME, utils.BASE, args.histo_path, args.base_name + "_histology"
    )
    regions = ["cervix", "cervical", "central", "ovarian"]
    horn = args.horn

    if not args.not_d:
        # If the dataset is downsampled
        uCT_directory = os.path.join(uCT_directory, "downsampled")

    # Add the muscle segmentation to the directory paths
    uCT_directory = os.path.join(uCT_directory, "muscle_segmentation")
    histo_directory = os.path.join(histo_directory, "muscle_segmentation")

    # Read data
    uCT_data = np.load(uCT_directory + "/angular_thickness.pkl",
                       allow_pickle=True)[horn]
    histo_data = np.load(histo_directory + "/angular_thickness.pkl",
                         allow_pickle=True)[horn]

    # Set nan values to 0
    uCT_nan_ind = np.where(np.isnan(uCT_data))
    histo_nan_ind = np.where(np.isnan(histo_data))
    uCT_data[uCT_nan_ind[0], uCT_nan_ind[1]] = 0
    histo_data[histo_nan_ind[0], histo_nan_ind[1]] = 0

    try:
        assert uCT_data.shape == histo_data.shape

    except AssertionError:
        sys.stderr.write(
            "Error: uCT and histology data show have same shape.\n")

    nb_samples = len(uCT_data)

    for i in range(uCT_data.shape[1]):
        t_stats = stats.ttest_rel(uCT_data[:, i], histo_data[:, i])
        p_stats = stats.pearsonr(uCT_data[:, i], histo_data[:, i])
        print(
            "{} section: r({}) = {:.2f}, p = {:.3f}".format(
                regions[i].capitalize(),
                nb_samples - 2,
                p_stats.statistic,
                p_stats.pvalue,
            )
        )
        print(
            "{} section: t({}) = {:.2f}, p = {:.3f}".format(
                regions[i].capitalize(),
                nb_samples - 2,
                t_stats.statistic,
                t_stats.pvalue,
            )
        )
