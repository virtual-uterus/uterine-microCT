#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
estrus-comp.py

Script to compare the results at different estrus stages
Author: Mathias Roesler
Date: 03/25
"""

import argparse
import os

import numpy as np
import thickness.utils as utils

from thickness.plots import plot_data
from thickness.constants import BASE, HOME

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compares the results between estrus stages"
    )
    parser.add_argument(
        "dir_path",
        type=str,
        metavar="dir-path",
        help="path from BASE to the dataset",
    )
    parser.add_argument(
        "estrus_config",
        type=str,
        metavar="estrus-config",
        help="name of the estrus configuration file",
    )
    parser.add_argument(
        "metric",
        type=str,
        choices=["muscle_thickness", "radius", "length"],
        help="name of the metric to use",
    )
    parser.add_argument(
        "--not-d",
        action="store_true",
        help="flag used if the dataset is not downsampled, default False",
    )

    # Parse input arguments
    args = parser.parse_args()

    load_directory = os.path.join(HOME, BASE, args.dir_path)
    param_file = os.path.join(load_directory, args.estrus_config + ".toml")

    # Load parameters
    params = utils.parse_TOML(param_file)
    datasets = params["phases"]  # Dataset names sorted by estrus
    metrics = dict()  # Empty dict to hold results

    for phase in datasets.keys():
        # Initalise metrics dict
        metrics[phase] = []

        for i, dataset in enumerate(datasets[phase]):
            # Create dataset specific variables
            base_name = dataset + "_PTA_1_Rec_Trans"
            data_directory = os.path.join(load_directory, base_name)

            if not args.not_d:
                # If the dataset is downsampled
                data_directory = os.path.join(data_directory, "downsampled")
                set_param_file = os.path.join(
                    data_directory, base_name + "_downsampled.toml"
                )

            else:
                # If not use top-level parameter file
                set_param_file = os.path.join(
                    data_directory,
                    base_name + ".toml",
                )

            set_params = utils.parse_TOML(set_param_file)
            split_nb = set_params["split_nb"]

            # Read metric data
            metric_directory = os.path.join(
                data_directory,
                "muscle_segmentation/",
            )
            metric_data = np.load(
                metric_directory + args.metric + ".pkl",
                allow_pickle=True,
            )

            if args.metric == "length":
                metrics[phase].append(
                    np.round(
                        [
                            np.mean(list(metric_data.values())),
                            np.std(list(metric_data.values())),
                        ],
                        2,
                    ),
                )
            else:
                mean_data = [
                    np.mean(list(metric_data.values())[0][split_nb:]),
                    np.mean(list(metric_data.values())[1][split_nb:]),
                ]
                std_data = [
                    np.std(list(metric_data.values())[0][split_nb:]),
                    np.std(list(metric_data.values())[1][split_nb:]),
                ]

                # Compute average std
                std_mean = np.sqrt(sum(np.power(std_data, 2)) / len(std_data))

                metrics[phase].append(
                    np.round([np.mean(mean_data), std_mean], 2),
                )

        metrics[phase] = np.array(metrics[phase])  # Convert to np array
    plot_data(metrics, args.metric)
