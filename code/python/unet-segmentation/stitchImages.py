#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# stitchImages.py: Script to stitch up images divided into blocks
# Author: Mathias Roesler
# Last modified: 04/24

import argparse
import glob
import sys
import os

import numpy as np
import skimage.io as skio
import utils.utils as utils

from skimage.util import img_as_ubyte

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Stiches back images that were divided into blocks"
    )

    parser.add_argument(
        "dir_path", type=str, metavar="dir-path",
        help="path from BASE to the dataset"
    )
    parser.add_argument(
        "base_name", type=str, metavar="base-name", help="name of the dataset"
    )
    parser.add_argument(
        "-e",
        "--extension",
        type=str,
        metavar="extension",
        help="extension for the images, default png",
        default="png",
    )
    parser.add_argument(
        "--not-d",
        action="store_true",
        help="flag used if the dataset is not downsampled, default False",
    )

    # Parse input arguments
    args = parser.parse_args()

    main_directory = os.path.join(utils.HOME, utils.BASE, args.dir_path,
                                  args.base_name)

    if not args.not_d:
        # If the dataset is downsampled
        main_directory = os.path.join(main_directory, "downsampled")
        param_file = os.path.join(
            main_directory, args.base_name + "_downsampled.toml")

    else:
        # If not use top-level parameter file
        param_file = os.path.join(main_directory, args.base_name + ".toml")

    load_directory = os.path.join(main_directory, "masks/")
    save_directory = os.path.join(main_directory, "stitched/")
    #
    # Check directory exists
    if not os.path.isdir(save_directory):
        os.mkdir(save_directory)

    # Get necessary parameters
    params = utils.parseTOML(param_file)
    original_width = params['nb_pixel_y']
    original_height = params['nb_pixel_x']
    image_name = params['prefix']

    img_list = sorted(glob.glob("*.{}".format(args.extension),
                                root_dir=load_directory))

    # Load first image to get height and width of the stack
    path = os.path.join(load_directory, img_list[0])

    if not os.path.isfile(path):
        sys.stderr.write("Error: {} is not a file.\n".format(path))
        exit()

    img = skio.imread(path, as_gray=True)
    height, width = img.shape

    # Get the padding values that were used for splitting
    pad_w = utils.findPadding(original_width, width)
    pad_h = utils.findPadding(original_height, height)

    # Calculate how many height and width blocks there are
    nb_h = (original_height + pad_h) // height
    nb_w = (original_width + pad_h) // width

    # Create temporary image for stitching
    tmp_img = np.zeros((original_height + pad_h, original_width + pad_w),
                       dtype=np.uint8)

    for i, img_name in enumerate(img_list):
        path = os.path.join(load_directory, img_name)

        if not os.path.isfile(path):
            sys.stderr.write("Error: {} is not a file.\n".format(path))
            exit()

        img = skio.imread(path, as_gray=True)

        # Find the limits of the block in the temporary image
        w_lims = [(i % nb_w) * width,
                  ((i % nb_w) + 1) * width]
        h_lims = [((i // nb_w) % nb_h) * height,
                  (((i // nb_w) % nb_h) + 1) * height]

        # Fill the correct block of the temporary image
        tmp_img[h_lims[0]:h_lims[1], w_lims[0]:w_lims[1]] = img

        # Save the image block
        if not (i+1) % (nb_h * nb_w):  # Account for 0 index by adding 1
            skio.imsave("{}/{}_{}.{}".format(
                save_directory, image_name,
                str("%03d" % ((i + 1) // (nb_h * nb_w))),
                args.extension),
                img_as_ubyte(tmp_img[:original_height,
                                     :original_width]),
                check_contrast=False)
