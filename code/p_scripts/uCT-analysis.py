#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# uCT-analysis.py: Script to analyse muscle thickness in uterine horns
# Author: Mathias Roesler
# Last modified: 06/23

import argparse
import os
import pickle

import numpy as np
import scipy.io

import thickness_analysis.plots as plots
import thickness_analysis.projection as projection
import utils.utils as utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Determines the thickness of the muscle layers "
        "from a uCT dataset"
    )

    parser.add_argument(
        "dir_path", type=str, metavar="dir-path", help="path from BASE "
        "to the dataset"
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

    # Parse input arguments
    args = parser.parse_args()

    load_directory = os.path.join(
        utils.HOME, utils.BASE, args.dir_path, args.base_name)

    if not args.not_d:
        # If the dataset is downsampled
        load_directory = os.path.join(load_directory, "downsampled")
        param_file = os.path.join(
            load_directory, args.base_name + "_downsampled.toml")

    else:
        # If not use top-level parameter file
        param_file = os.path.join(load_directory, args.base_name + ".toml")

    # Load parameters
    params = utils.parseTOML(param_file)
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
    avg_slice_thickness = dict()
    errors = dict()

    for i, horn in enumerate(horns):
        print("Processing {} horn".format(horn))
        print("   Loading mask stack")
        mask_stack = utils.loadImageStack(
            os.path.join(
                load_directory, "{}".format(horn)), extension=args.extension
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
        muscle_thickness, slice_thickness, radius = \
            projection.estimateMuscleThickness(
                mask_stack, centreline, args.points,
                params[horn]["slice_nbs"], horn
            )

        # Estimate horn length
        print("   Estimating horn length")
        centreline_dict = scipy.io.loadmat(
            load_directory + "/centreline.mat"
        )
        centreline = np.transpose(centreline_dict["centreline"])
        centreline = np.round(centreline).astype(int)  # Convert to int

        match horn:
            case "left":
                ind = np.where(centreline[:nb_slices, 0:4] == 0)[0]
                horn_start = np.append(centreline[ind[0], 0:2], 0)
                # Divide len by 2 because np.where doubles the length
                horn_end = np.append(centreline[ind[-1], 0:2], len(ind) / 2)

            case "right":
                ind = np.where(centreline[:nb_slices, 2:6] == 0)[0]
                horn_start = np.append(centreline[ind[0], 4:6], 0)
                # Divide len by 2 because np.where doubles the length
                horn_end = np.append(centreline[ind[-1], 4:6], len(ind) / 2)

        length = np.linalg.norm(horn_end - horn_start)

        # Rescale the thickness to mm
        muscle_thickness *= params["scaling_factor"]
        slice_thickness *= params["scaling_factor"]
        radius *= params["scaling_factor"]
        length *= params["scaling_factor"]

        print(
            "{} horn muscle thickness: {:.2f} \u00B1 {:.2f} mm".format(
                horn, np.mean(muscle_thickness), np.std(muscle_thickness)
            )
        )
        print(
            "{} horn radius: {:.2f} \u00B1 {:.2f} mm".format(
                horn, np.mean(radius), np.std(radius)
            )
        )
        print(
            "{} horn length: {:.2f} mm".format(
                horn, length
            )
        )

        if args.switch:
            avg_thickness[horns[i - 1]] = utils.movingAverage(
                muscle_thickness, muscle_win_size
            ).round(5)
            avg_slice_thickness[horns[i - 1]] = utils.circularAverage(
                slice_thickness, circular_win_size
            ).round(5)
            errors[horns[i - 1]] = utils.movingStd(muscle_thickness,
                                                   std_win_size)

        else:
            avg_thickness[horn] = utils.movingAverage(
                muscle_thickness, muscle_win_size).round(5)
            avg_slice_thickness[horn] = utils.circularAverage(
                slice_thickness, circular_win_size).round(5)
            errors[horn] = utils.movingStd(muscle_thickness, std_win_size)

    # Save angular thickness
    with open(load_directory + "/angular_thickness.pkl", "wb") as f:
        pickle.dump(avg_slice_thickness, f)

    # Save muscle thickness
    with open(load_directory + "/muscle_thickness.pkl", "wb") as f:
        pickle.dump(avg_thickness, f)

    # Plot everything
    plots.plotMuscleThickness(avg_thickness, errors)

    if len(horns) == 2:
        for horn in horns:
            plots.plotAngularThickness(
                {horn: avg_slice_thickness[horn]}, projection=args.polar
            )
