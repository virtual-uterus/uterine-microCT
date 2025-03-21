#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
uCT-analysis.py

Script to analyse muscle thickness in uterine horns
Author: Mathias Roesler
Date: 06/23
"""

import argparse
import os
import pickle

import numpy as np
import scipy.io

import thickness.plots as plots
import thickness.projection as projection
import thickness.utils as utils

from thickness.constants import BASE, HOME

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Determines the thickness of the myometrium from uCT data"
    )

    parser.add_argument(
        "dir_path",
        type=str,
        metavar="dir-path",
        help="path from BASE to the dataset",
    )
    parser.add_argument(
        "base_name", type=str, metavar="base-name", help="name of the dataset"
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
        choices={"left", "right", "both"},
        help="horn to process",
        default="both",
    )
    parser.add_argument(
        "-p",
        "--points",
        type=int,
        metavar="points",
        help="number of points to use for the projection, default 128",
        default=128,
    )
    parser.add_argument(
        "-P",
        "--polar",
        action="store_true",
        help="flag used to plot the angular thickness in polar projection, "
        "default False",
    )
    parser.add_argument(
        "-s",
        "--switch",
        action="store_true",
        help="switches the labels of the left and right horn, default False",
    )
    parser.add_argument(
        "--not-d",
        action="store_true",
        help="flag used if the dataset is not downsampled, default False",
    )
    parser.add_argument(
        "--no-plot",
        action="store_false",
        help="flag used to plot the results, default False",
    )

    # Parse input arguments
    args = parser.parse_args()

    load_directory = os.path.join(
        HOME,
        BASE,
        args.dir_path,
        args.base_name,
    )

    if not args.not_d:
        # If the dataset is downsampled
        load_directory = os.path.join(load_directory, "downsampled")
        param_file = os.path.join(
            load_directory,
            args.base_name + "_downsampled.toml",
        )

    else:
        # If not use top-level parameter file
        param_file = os.path.join(load_directory, args.base_name + ".toml")

    # Load parameters
    params = utils.parse_TOML(param_file)
    split_nb = params["split_nb"]  # Get horn separation slice
    weight = params["weight"] * 1e-3  # Weight in mg for normalisation
    params = params["thickness"]  # Extract the thickness parameters

    # Add the muscle segmentation to the load directory
    load_directory = os.path.join(load_directory, "muscle_segmentation")

    # Convert both to left and right
    if args.horn == "both":
        horns = ["left", "right"]

    else:
        horns = [args.horn]

    # Dicts for results of both horns
    avg_thickness = dict()
    normalised_thickness = dict()
    avg_slice_thickness = dict()
    errors = dict()
    radius_dict = dict()
    length_dict = dict()

    for i, horn in enumerate(horns):
        if args.switch:
            print_horn = horns[i - 1]

        else:
            print_horn = horn

        print("Processing {} horn".format(print_horn))
        print("   Loading mask stack")
        mask_stack = utils.load_image_stack(
            os.path.join(load_directory, "{}".format(horn)),
            extension=args.extension,
        )

        nb_imgs = len(mask_stack)

        # Window sizes for different moving averages
        muscle_win_size = round(0.10 * nb_imgs)
        std_win_size = 20
        circular_win_size = round(0.04 * args.points)

        print("   Loading centreline")
        centreline_dict = scipy.io.loadmat(
            load_directory + "/{}/centreline.mat".format(horn)
        )
        centreline = np.transpose(centreline_dict["centreline"])
        centreline = np.round(centreline).astype(int)  # Convert to int
        nb_slices = len(centreline) - 1  # Number of slices in the horn

        print("   Estimating muscle thickness")
        muscle_thickness, slice_thickness, radius = (
            projection.estimate_muscle_thickness(
                mask_stack,
                centreline,
                args.points,
                params[horn]["slice_nbs"],
                horn,
            )
        )

        # Estimate horn length
        print("   Estimating horn length")
        centreline_dict = scipy.io.loadmat(load_directory + "/centreline.mat")
        centreline = np.transpose(centreline_dict["centreline"])

        match horn:
            case "left":
                ind = np.where(centreline[:nb_slices, 0] != 0)[0]
                last_slice = ind[-1]  # Get the last slice index of horn
                centre_vectors = np.diff(centreline[:last_slice, 0:2], axis=0)
            case "right":
                ind = np.where(centreline[:nb_slices, 5] != 0)[0]
                last_slice = ind[-1]  # Get the last slice index of horn
                centre_vectors = np.diff(centreline[:last_slice, 4:6], axis=0)

        coordinates = np.ones(
            (last_slice - split_nb - 1, 3),
        )  # Account for diff
        # Add the centre vector coordinates
        coordinates[:, 0:2] = centre_vectors[split_nb:]
        length = np.sum(np.linalg.norm(coordinates, axis=1))

        # Rescale the thickness to mm
        muscle_thickness *= params["scaling_factor"]
        slice_thickness *= params["scaling_factor"]
        radius *= params["scaling_factor"]
        length *= params["scaling_factor"]

        # Populate dictionnaries for pickling data
        avg_thickness[print_horn] = utils.moving_average(
            muscle_thickness, muscle_win_size
        ).round(5)
        avg_slice_thickness[print_horn] = utils.circular_average(
            slice_thickness, circular_win_size
        ).round(5)
        errors[print_horn] = utils.moving_std(muscle_thickness, std_win_size)
        radius_dict[print_horn] = radius
        length_dict[print_horn] = length

        print(
            "{} horn muscle thickness: {:.2f} \u00b1 {:.2f} mm".format(
                print_horn,
                np.mean(muscle_thickness[split_nb:]),
                np.std(muscle_thickness[split_nb:]),
            )
        )
        print(
            "{} horn radius: {:.2f} \u00b1 {:.2f} mm".format(
                print_horn,
                np.mean(radius[split_nb:]),
                np.std(radius[split_nb:]),
            )
        )
        print("{} horn length: {:.2f} mm".format(print_horn, length))

        # Normalise by weight
        normalised_thickness[print_horn] = avg_thickness[print_horn] / weight
        radius_dict[print_horn] /= weight

    # Save angular thickness
    with open(load_directory + "/angular_thickness.pkl", "wb") as f:
        pickle.dump(avg_slice_thickness, f)

    # Save muscle thickness
    with open(load_directory + "/muscle_thickness.pkl", "wb") as f:
        pickle.dump(normalised_thickness, f)

    # Save radius
    with open(load_directory + "/radius.pkl", "wb") as f:
        pickle.dump(radius_dict, f)

    # Save horn length
    with open(load_directory + "/length.pkl", "wb") as f:
        pickle.dump(length_dict, f)

    # Plot everything
    if args.no_plot:
        plots.plot_muscle_thickness(avg_thickness, errors)

        if len(horns) == 2:
            for horn in horns:
                plots.plot_angular_thickness(
                    {horn: avg_slice_thickness[horn]}, projection=args.polar
                )
