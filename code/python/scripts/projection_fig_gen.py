#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# projection_fig_gen.py: Script to generate the projection point figure
# Author: Mathias Roesler
# Last modified: 06/23

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.io

import thickness_analysis.plots as plots
import thickness_analysis.projection as projection
import utils.utils as utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Specific script to generate the projection point plot"
    )

    parser.add_argument(
        "dir_path", type=str, metavar="dir-path",
        help="path from BASE to the dataset"
    )
    parser.add_argument(
        "base_name", type=str, metavar="base-name", help="name of the dataset"
    )
    parser.add_argument(
        "img_nb", type=str, metavar="img-nb", help="number of the image to use"
    )
    parser.add_argument(
        "-e",
        "--extension",
        type=str,
        metavar="extension",
        help="extension for the saved images, default png",
        default="png",
    )
    parser.add_argument(
        "--horn",
        type=str,
        choices={"left", "right"},
        help="horn to process, default right",
        default="right",
    )
    parser.add_argument(
        "-p",
        "--points",
        type=int,
        metavar="points",
        help="number of points to use for the projection, default 8",
        default=8,
    )
    parser.add_argument(
        "--not-d",
        action="store_true",
        help="flag used if the dataset is not downsampled, default False",
    )

    # Parse input arguments
    args = parser.parse_args()

    load_directory = os.path.join(utils.HOME, utils.BASE, args.dir_path,
                                  args.base_name)

    if not args.not_d:
        # If the dataset is downsampled
        load_directory = os.path.join(load_directory, "downsampled")
        param_file = os.path.join(load_directory,
                                  args.base_name + "_downsampled.toml")

    else:
        # If not use top-level parameter file
        param_file = os.path.join(load_directory, args.base_name + ".toml")

    # Load parameters
    params = utils.parseTOML(param_file)

    # Get the original image
    original_img_name = os.path.join(
        load_directory,
        params["prefix"] + "_" + args.img_nb + "." + args.extension
    )

    # Add the muscle segmentation to the load directory
    load_directory = os.path.join(load_directory, "muscle_segmentation")

    # Get the original mask
    original_mask_name = os.path.join(
        load_directory,
        params["prefix"] + "_" + args.img_nb + "." + args.extension
    )

    load_directory = os.path.join(load_directory, args.horn)

    # Load the centreline
    centreline_dict = scipy.io.loadmat(load_directory + "/centreline.mat")
    centreline = np.transpose(centreline_dict["centreline"])
    centreline = np.round(centreline).astype(int)  # Round and convert to int

    # Image to use for projection
    rotated_mask_name = os.path.join(
        load_directory,
        params["prefix"] + "_" + args.img_nb + "." + args.extension
    )

    # Load all the images
    original_img = plt.imread(original_img_name)
    original_mask = plt.imread(original_mask_name)
    rotated_mask = plt.imread(rotated_mask_name)

    # Rectify image number because index starts at 0
    projection_points = projection.findProjectionPoints(
        rotated_mask, centreline[int(args.img_nb) - 1], args.points, args.horn
    )

    # Plot everything
    fig, ax = plt.subplots(dpi=300)
    plt.imshow(original_img, cmap="gray")
    fig, ax = plt.subplots(dpi=300)
    plt.imshow(original_mask, cmap="gray")
    fig, ax = plt.subplots(dpi=300)
    plt.imshow(rotated_mask, cmap="gray")
    plots.plotProjectionPoints(
        rotated_mask, centreline[int(args.img_nb) - 1], projection_points
    )
