#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nifti-converter.py

Programs that converts a stack of images into a
nifti volume. Requires SimpleITK package.
Author: Mathias Roesler
Date: 06/23
"""

import argparse
import glob
import os
import sys

import SimpleITK as sitk

from thickness.constants import BASE, HOME

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Converts uCT dataset into the nifti format"
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
        help="extension of the images, default png",
        default="png",
    )
    parser.add_argument(
        "--segmentation",
        type=str,
        default="",
        help="name of the segmentation folder, default ''",
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

    # Add the segmentation folder if provided
    load_directory = os.path.join(load_directory, args.segmentation)

    # Get the paths of the images in the load directory
    img_path = os.path.join(load_directory, "*." + args.extension)

    if not os.path.exists(load_directory):
        sys.stderr.write(
            f"Error: the input folder {load_directory} does not exist\n",
        )
        exit()

    # Sort the images
    img_list = sorted(glob.glob(img_path))

    if len(img_list) == 0:
        sys.stderr.write(
            f"Error: the folder {load_directory} does not contain any images\n"
        )
        exit()

    # Add the .nii.gz extension to the base name
    if len(args.segmentation) > 0:
        save_name = args.base_name + "_" + args.segmentation + ".nii.gz"

    else:
        save_name = args.base_name + ".nii.gz"

    # Read in all the images
    reader = sitk.ImageSeriesReader()
    reader.SetFileNames(img_list)
    volume = reader.Execute()

    sitk.WriteImage(volume, os.path.join(load_directory, save_name))
