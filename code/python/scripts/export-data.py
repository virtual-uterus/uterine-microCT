#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
export-data.py

Exports data to R for statistical analysis
Author: Mathias Roesler
Date: 09/25
"""

import os
import sys
import argparse

import numpy as np
import pandas as pd
import thickness.utils as utils

from thickness.constants import BASE, HOME

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Export data to R for statistical analysis"
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
        choices=["muscle_thickness", "radius", "length", "endometrium_volume"],
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

    phase_list = []  # Empty list to contain the phases
    dataset_list = []  # Empty list to contain dataset ids
    value_list = []  # Empty list to contain average values

    for phase in datasets.keys():
        for i, dataset in enumerate(datasets[phase]):
            phase_list.append(phase)  # Phase list update
            dataset_list.append(dataset)  # Add datasets in the correct order

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
            if args.metric == "endometrium_volume":
                metric_directory = os.path.join(
                    data_directory,
                    "endometrium_segmentation/",
                )
            else:
                metric_directory = os.path.join(
                    data_directory,
                    "muscle_segmentation/",
                )
            metric_data = np.load(
                metric_directory + args.metric + ".pkl",
                allow_pickle=True,
            )

            if args.metric == "length" or args.metric == "endometrium_volume":
                value_list.append(np.mean(list(metric_data.values())))
            else:
                value_list.append(
                    np.mean(
                        [
                            np.mean(list(metric_data.values())[0][split_nb:]),
                            np.mean(list(metric_data.values())[1][split_nb:]),
                        ]
                    )
                )

    # Create data frame and export to csv for R
    df = pd.DataFrame(
        {"Phase": phase_list, "Experiment": dataset_list, "Value": value_list},
    )

    try:
        df.to_csv(
            load_directory + "/exports/" + args.metric + ".csv",
            index=False,
        )

    except OSError as e:
        sys.stderr.write("Error: {}".format(e))
        exit()
