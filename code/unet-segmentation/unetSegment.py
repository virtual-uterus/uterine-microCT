#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# unetSegment.py: Script to segment microCT using unet model
# Author: Mathias Roesler
# Last modified: 05/24

import os
import glob
import argparse

import numpy as np
import tensorflow as tf
import utils.utils as utils
import skimage.io as skio

from skimage.filters import threshold_otsu
from skimage.util import img_as_ubyte

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Segment microCT dataset with unet model"
    )

    parser.add_argument(
        "dir_path", type=str, metavar="dir-path",
        help="path from BASE to the dataset"
    )
    parser.add_argument(
        "base_name", type=str, metavar="base-name", help="name of the dataset"
    )
    parser.add_argument(
        "model", type=str, metavar="model",
        help="Saved model to use"
    )
    parser.add_argument(
        "--not-d",
        action="store_true",
        help="flag used if the dataset is not downsampled, default False",
    )
    parser.add_argument(
        "-e",
        "--extension",
        type=str,
        metavar="extension",
        help="extension for the saved images, default png",
        default="png",
    )

    # Parse input arguments
    args = parser.parse_args()

    main_directory = os.path.join(utils.HOME, utils.BASE, args.dir_path,
                                  args.base_name)

    if not args.not_d:
        # If the dataset is downsampled
        main_directory = os.path.join(main_directory, "downsampled")

    save_directory = os.path.join(main_directory, "masks")
    load_directory = os.path.join(main_directory, "imgs")

    # Check directory exists
    if not os.path.isdir(save_directory):
        os.mkdir(save_directory)

    imgs = utils.loadImageStack(load_directory)
    img_names = sorted(glob.glob("*.{}".format(args.extension),
                                 root_dir=load_directory))

    # Convert to floats between 0 and 1
    imgs = np.asarray(imgs, dtype=np.float32) / imgs.max()

    # Ensure shape are correct
    if len(imgs.shape) < 4:
        imgs = imgs.reshape(imgs.shape[0],
                            imgs.shape[1],
                            imgs.shape[2],
                            1)

    # Load trained model
    model = tf.keras.models.load_model(args.model)
    masks = model.predict(imgs)

    # Reset max value to 255
    masks *= 255
    # Convert masks to uint8
    masks = masks.astype('uint8')

    for i in range(masks.shape[0]):
        # Get mask
        mask = np.squeeze(masks[i, :, :])

        # Binarise image
        threshold = threshold_otsu(mask)
        mask = mask > threshold

        skio.imsave("{}/{}.{}".format(
            save_directory, os.path.splitext(img_names[i])[0],
            args.extension), img_as_ubyte(mask), check_contrast=False)
