#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# mesh-converter.py: Script that converts a vtu mesh to cmgui format
# Author: Mathias Roesler
# Last modified: 06/23

import argparse
import os
import sys

import meshio

import utils.utils as utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Converts an annotated vtu or vtk mesh to a cmgui format"
    )

    # Parse input arguments
    parser.add_argument(
        "mesh_name", type=str, metavar="mesh-name",
        help="name of the mesh to convert"
    )
    parser.add_argument(
        "--mesh-dir",
        type=str,
        default="mesh/",
        help="path from BASE to the mesh, default mesh/",
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
    mesh_path = os.path.join(utils.HOME, utils.BASE, args.mesh_dir)
    mesh_file = mesh_path + "/" + args.mesh_name

    # Read the mesh file
    print("Loading mesh {}".format(mesh_file + "." + args.extension))
    mesh = meshio.read(mesh_file + "." + args.extension)

    # Extract information
    nodes = mesh.points

    try:
        thickness = mesh.point_data["thickness"]
        thickness_flag = True

    except KeyError:
        sys.stderr.write("Warning: no thickness data found\n")
        thickness = None
        thickness_flag = False

    if len(thickness.shape) == 2:
        # Using a vtu format that needs to be reshaped
        thickness = thickness.reshape(thickness.shape[0])

    try:
        elements = mesh.cells_dict["tetra"]
        vol_flag = True

    except KeyError:
        elements = mesh.cells_dict["triangle"]
        vol_flag = False

    # Write EX files
    print("Writing exnode file")
    utils.writeExNode(mesh_file + ".exnode", nodes, thickness)

    print("Writing exelem file")
    if vol_flag:
        utils.writeExElemVol(mesh_file + ".exelem", elements, thickness_flag)

    else:
        utils.writeExElemSurf(mesh_file + ".exelem", elements, thickness_flag)
