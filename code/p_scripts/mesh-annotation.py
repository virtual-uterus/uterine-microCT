#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# mesh-annotation.py: Script that adds the thickness data to a vtu/vtk mesh
# Author: Mathias Roesler
# Last modified: 06/23

import argparse
import os

import meshio
import numpy as np
import scipy.io

import utils.utils as utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotates a vtu or vtk mesh with the thickness values "
        "of each slice"
    )

    # Parse input arguments
    parser.add_argument(
        "base_name", type=str, metavar="base-name", help="name of the dataset"
    )
    parser.add_argument(
        "--mesh-dir",
        type=str,
        default="mesh/",
        help="path from BASE to the mesh, default mesh/",
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="microCT/data",
        help="path form BASE to the thickness data, default microCT/data",
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
    parser.add_argument(
        "-e",
        "--extension",
        choices={"vtu", "vtk"},
        help="mesh extesion, default value vtk",
        default="vtk",
    )

    args = parser.parse_args()

    # Set arguments
    mesh_directory = os.path.join(utils.HOME, utils.BASE, args.mesh_dir)
    thickness_directory = os.path.join(
        utils.HOME, utils.BASE, args.data_dir, args.base_name
    )
    mesh_name = os.path.join(
        mesh_directory, args.base_name + "_volumetric_mesh")

    if not args.not_d:
        # If the dataset is downsampled
        thickness_directory = os.path.join(thickness_directory, "downsampled")
        param_file = os.path.join(
            thickness_directory, args.base_name + "_downsampled.toml"
        )

    else:
        param_file = os.path.join(
            thickness_directory, args.base_name + ".toml")

    # Load parameters
    params = utils.parseTOML(param_file)
    nb_used_slices = params["nb_used_slices"]  # Get number of slices to use
    start_nb = params["thickness"]["start_nb"]  # Number of slices not rotated

    # Add the muscle segmentation to the load directory
    thickness_directory = os.path.join(
        thickness_directory, "muscle_segmentation")

    # Convert both to left and right
    if args.horn == "both":
        horns = ["left", "right"]

    else:
        horns = [args.horn]

    # Get the centreline
    centreline_dict = scipy.io.loadmat(thickness_directory + "/centreline.mat")
    centreline = np.transpose(centreline_dict["centreline"])
    centreline = np.round(centreline).astype(int)  # Round and convert to int

    # Read the mesh file
    print("Loading mesh {}".format(mesh_name + "." + args.extension))
    mesh = meshio.read(mesh_name + "." + args.extension)
    nb_points = len(mesh.points)

    # Re-centre z coordinates of the mesh with origin
    mesh.points[:, 2] = mesh.points[:, 2] - min(mesh.points[:, 2])
    mesh.points[:, 0] = -mesh.points[:, 0]
    mesh.points[:, 1] = -mesh.points[:, 1]

    # Read thickness data
    thickness_data = np.load(
        thickness_directory + "/muscle_thickness.pkl", allow_pickle=True
    )

    # Thickness data dictionary to be added to the mesh
    point_data_dict = dict()
    point_data_array = np.zeros((nb_points, 1))
    point_data_name = "thickness"

    # Split the mesh elements into left and right
    element_idx_dict = dict()
    x_delta = mesh.points[:, 0].max() - mesh.points[:, 0].min()

    if len(horns) == 1:
        # If single horn all elements are associated with horn
        element_idx_dict[horns[0]] = np.arange(len(mesh.points))

    else:
        # If two horns split elements into left and right
        element_idx_dict["left"] = np.where(
            abs(mesh.points[:, 0]) < x_delta)[0]
        element_idx_dict["right"] = np.where(
            abs(mesh.points[:, 0]) >= x_delta)[0]

    # Split the centre points into left and right
    centrepoint_coords = dict()
    centrepoint_coords["left"] = centreline[:, 0:2]
    centrepoint_coords["right"] = centreline[:, 4:6]

    for i, horn in enumerate(horns):
        print("Annotating {} horn".format(horn))
        if args.switch:
            # If the centrepoints need to be switched
            i = i - 1

        thickness = thickness_data[horn]
        element_idx = element_idx_dict[horns[i]]
        elements = mesh.points[element_idx]
        centrepoints = centrepoint_coords[horns[i]]

        for k in range(len(thickness)):
            # Find the centre vector
            if k <= start_nb:
                # Unrotated slices
                centre_vector = np.array([0, 0, 1])

            else:
                if k >= len(thickness) - nb_used_slices:
                    # Use less slices to get centre vector
                    next_centrepoint = centrepoints[len(thickness) - 1]
                    z_centre = len(thickness) - k
                else:
                    next_centrepoint = centrepoints[k + nb_used_slices]
                    z_centre = nb_used_slices

                centre_vector = next_centrepoint - centrepoints[k]

                # Add z component
                centre_vector = np.append(centre_vector, z_centre)

                # Ensure that the direction of centre vector is consistent
                condition_1 = centre_vector[0] < 0 and horns[i] == "right"
                condition_2 = centre_vector[0] > 0 and horns[i] == "left"

                if condition_1 or condition_2:
                    centre_vector[0] = -centre_vector[0]

            # Slice plane normal vector and origin point
            plane_normal = centre_vector / np.linalg.norm(centre_vector)
            plane_origin = np.array([centrepoints[k, 0], centrepoints[k, 1],
                                     k])

            norms = np.dot(elements - plane_origin, plane_normal)
            idx_list = element_idx[(norms <= 15) & (norms >= 0)]
            point_data_array[idx_list] = round(thickness[k], 5)

    # Add the data dictionary to the mesh
    point_data_dict[point_data_name] = point_data_array
    mesh.point_data = point_data_dict

    # Flip mesh back to original position
    mesh.points[:, 0] = -mesh.points[:, 0]
    mesh.points[:, 1] = -mesh.points[:, 1]

    # Save new mesh
    print("Saving mesh {}".format(mesh_name + "_annotated." + args.extension))
    mesh.write(mesh_name + "_annotated." + args.extension)
