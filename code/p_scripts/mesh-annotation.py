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


def getIndices(normal_vector, elements, plane_distance, centre_norm):
    """Finds the indices of the elements to annotate with the thickness.

    Arguments:
    normal_vector -- ndarray, normal vector for the plane.
    elements -- ndarray, list of elements in the mesh to sort.
    plane_distance -- int, distance between origin and cut planes.
    centre_norm -- float, norm of the normal vector.

    Return:
    idx_list -- list[float], list of indices of the correct elements.

    """
    ele_dot_prod = np.dot(normal_vector, elements)

    # Get points between the two planes
    dist_to_first_plane = (ele_dot_prod + plane_distance) / centre_norm
    dist_to_second_plane = (ele_dot_prod - plane_distance) / centre_norm

    idx_list = [
        idx
        for idx, (a, b) in enumerate(
            zip(
                dist_to_first_plane > 0.0,
                dist_to_second_plane < 0.0,
            )
        )
        if a and b
    ]

    return idx_list


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

    tolerance = 5

    for i, horn in enumerate(horns):
        print("Annotating {} horn".format(horn))
        if args.switch:
            # If the centrepoints need to be switched
            i = i - 1

        thickness = thickness_data[horn]

        if horns[i] == "left":
            centrepoints = centreline[0: len(thickness), 0:2]
        else:
            centrepoints = centreline[0: len(thickness), 4:6]

        centrepoints_diff = np.diff(centrepoints, axis=0)
        z_components = np.ones((len(centrepoints_diff), 1))

        # Create normalised centre vectors
        centre_vectors = np.append(centrepoints_diff, z_components, axis=1)
        centre_norms = np.repeat(
            np.linalg.norm(centre_vectors, axis=1), 3
        )  # Repeat to be able to reshape
        centre_norms = np.reshape(
            centre_norms, centre_vectors.shape
        )  # Reshape for division
        centre_vectors_norm = centre_vectors / centre_norms

        for j in range(len(thickness)):
            x_split = (centreline[j, 0] + centreline[j, 4]) / 2

            if j >= len(centre_vectors_norm):
                centre_vector = np.array([0, 0, 1])
                centre_norm = 1
            else:
                centre_vector = centre_vectors_norm[j]
                centre_norm = centre_norms[j][0]

            # Get the x limited indices
            if horns[i] == "left":
                x_idx = np.where(mesh.points[:, 0] < x_split)[0]
            else:
                x_idx = np.where(mesh.points[:, 0] >= x_split)[0]

            elements = mesh.points[x_idx] - np.append(
                centrepoints[j], j
            )  # Reduced set of points and recentre

            idx_list = getIndices(
                centre_vector, np.transpose(elements), plane_distance, centre_norm
            )

            if j <= 20 or j >= len(thickness) - 20:
                extra_idx = getIndices(
                    np.array([0, 0, 1]), np.transpose(elements), plane_distance, 1
                )
                idx_list = np.append(idx_list, extra_idx)
            point_data_array[x_idx[idx_list]] = round(thickness[j], 3)

    # Add the data dictionary to the mesh
    point_data_dict[point_data_name] = point_data_array
    mesh.point_data = point_data_dict

    # Flip mesh back to original position
    mesh.points[:, 0] = -mesh.points[:, 0]
    mesh.points[:, 1] = -mesh.points[:, 1]

    # Save new mesh
    print("Saving mesh {}".format(mesh_name + "_annotated." + args.extension))
    mesh.write(mesh_name + "_annotated." + args.extension)
