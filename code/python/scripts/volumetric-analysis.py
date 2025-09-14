#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
volumetric-analysis.py

Script to analyse endometrium and myometrium volume in uterine horns
Author: Mathias Roesler
Date: 07/25
"""

import argparse
import os
import pickle

import numpy as np

import thickness.utils as utils

from thickness.constants import BASE, HOME

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Determines endometrium and myometrium volume from uCT"
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
    resolution = params["resolution"] * 1e-3  # In mm

    data_type = ["endometrium", "muscle"]

    for data in data_type:
        # Add the data_type segmentation to the load directory
        data_load_directory = os.path.join(
            load_directory,
            data + "_segmentation",
        )

        # Convert both to left and right
        if args.horn == "both":
            horns = ["left", "right"]

        else:
            horns = [args.horn]

        volume_dict = dict()

        for i, horn in enumerate(horns):
            if args.switch:
                print_horn = horns[i - 1]

            else:
                print_horn = horn

            print("Processing {} horn".format(print_horn))
            print("   Loading mask stack")
            mask_stack = utils.load_image_stack(
                os.path.join(data_load_directory, "{}".format(horn)),
                extension=args.extension,
            )

            nb_pixels = np.sum(mask_stack[:, :, split_nb:] >= 1)
            volume = nb_pixels * (resolution**3)
            volume_dict[print_horn] = volume / weight

            print(
                "{} horn {} volume: {:.2f} mm3/mg".format(
                    print_horn,
                    data,
                    volume / weight,
                )
            )

        with open(data_load_directory + "/" + data + "_volume.pkl", "wb") as f:
            pickle.dump(volume_dict, f)
