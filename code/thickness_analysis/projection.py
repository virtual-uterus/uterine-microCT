#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# projection.py: Functions used to find the projections points
# Author: Mathias Roesler
# Last modified: 04/24

import logging
import sys

import numpy as np
import skimage.draw as skd
import utils.utils as utils


def findLineCoordinates(img_shape, centre_point, theta):
    """Finds the line in the image given an angle and a point.

    Arguments:
    img_shape -- ndarray, shape of the image that is being analysed.
    centre_point -- ndarray, coordinates of the centre point (XY).
    theta -- angle in radians.

    Return:
    line_x -- list[int], x coordinates of the line that belong to the
        image.
    line_y -- list[int], y coordinates of the line that belong to the
        image.

    """
    # Ensure that coordinates are integers
    x_centre = int(np.round(centre_point[0]))
    y_centre = int(np.round(centre_point[1]))

    if theta != np.pi:
        y_points = np.arange(0, img_shape[0], 0.3)
        x_points = ((y_centre - y_points) * np.cos(theta)) / np.sin(
            theta) + x_centre
        x_points = x_points.astype(int)  # Convert to int
        y_points = y_points.astype(int)  # Convert to int

        # Find the points that are in the image
        intersection = np.intersect1d(
            np.where(x_points >= 0), np.where(x_points < img_shape[1])
        )
        x_points = x_points[intersection]
        y_points = y_points[intersection]

    else:
        # Case of the vertical line that goes through the centre
        x_points = np.arange(img_shape[1], dtype=int)
        y_points = np.ones(img_shape[1], dtype=int) * y_centre

    # Create the line and get the pixel values that it cuts through
    line_y, line_x = skd.line(y_points[0], x_points[0], y_points[-1],
                              x_points[-1])

    return line_x, line_y


def separationLine(img_shape, centre_point, normal):
    """Finds the line in the image given a normal and a point.

    Arguments:
    img_shape -- ndarray, shape of the image that is being analysed.
    centre_point -- ndarray, coordinates of the centre point (XY).
    normal -- ndarray, coordinates of the normal vector.

    Return:
    line_x -- list[int], x coordinates of the line that belong to the
        image.
    line_y -- list[int], y coordinates of the line that belong to the
        image.

    """
    # Ensure that coordinates are integers
    x_centre = int(np.round(centre_point[0]))
    y_centre = int(np.round(centre_point[1]))

    t = np.arange(-max(img_shape), max(img_shape))
    x_points = x_centre + normal[0] * t
    y_points = y_centre + normal[1] * t

    # Create the line and get the pixel values that it cuts through
    points = [
        (x, y)
        for x, y in zip(x_points, y_points)
        if 0 <= x < img_shape[1] and 0 <= y < img_shape[0]
    ]
    x_points, y_points = zip(*points)
    line_y, line_x = skd.line(
        int(y_points[0]), int(x_points[0]), int(y_points[-1]),
        int(x_points[-1])
    )

    return line_x, line_y


def excludeCentralPoints(img, centre_points, projection_points, horn):
    """Changes coordinates of points in the exclusion zone to [0, 0]

    Arguments:
    img -- ndarray, image to analyse.
    centre_points -- list[int], coordinates of the centre points (XY).
    projection_points -- ndarray, list of the coordinates of the
        projection points.
    horn -- str {left, right}, horn that is being analysed.

    Return:
    projection_points -- ndarray, list of the coordinates of the
        projection points with [0, 0] coordinates if in the exclusion zone.
    """
    # Ensure that coordinates are integers
    centre_points = np.array([np.round(x) for x in centre_points], dtype=int)

    # Create the line between left and right horn centre points
    hline_y, hline_x = skd.line(centre_points[1], centre_points[0],
                                centre_points[5], centre_points[4])
    hline = img[hline_y, hline_x]

    # Find the indices of the inner and outer edges of the cervix
    coords = np.where(hline[:-1] != hline[1:])[0]

    # Normal vector for the horizontal line
    normal = utils.getVector(
        np.array([centre_points[5], centre_points[0]]),
        np.array([centre_points[1], centre_points[4]]),
    )

    match horn:
        case "left":
            h_point = [hline_x[coords[0]-1], hline_y[coords[0]-1]]

            # Create the lines delimiting the exclusion zone
            first_vline_x, first_vline_y = separationLine(
                img.shape, h_point, normal)
            second_vline_x, second_vline_y = separationLine(
                img.shape, centre_points[2:4], normal)

        case "right":
            h_point = [hline_x[coords[1]+1], hline_y[coords[1]+1]]

            # Create the lines delimiting the exclusion zone
            first_vline_x, first_vline_y = separationLine(
                img.shape, centre_points[2:4], normal)
            second_vline_x, second_vline_y = separationLine(
                img.shape, h_point, normal)

    first_vline_vector = utils.getVector(
        np.array([first_vline_x[0], first_vline_y[0]]),
        np.array([first_vline_x[-1], first_vline_y[-1]]))
    second_vline_vector = utils.getVector(
        np.array([second_vline_x[0], second_vline_y[0]]),
        np.array([second_vline_x[-1], second_vline_y[-1]]))

    unit_vector = np.array([1, 0])  # Unit vector x direction

    # Get angles between x direction and exclusion lines
    first_angle = utils.getAngle(first_vline_vector, unit_vector)
    second_angle = utils.getAngle(second_vline_vector, unit_vector)

    # Create an array to store point indices
    ind_to_del = np.array([], dtype=int)

    for i, point in enumerate(projection_points):
        # Check if points are in exclusion zone
        vec1 = utils.getVector(
            np.array([first_vline_x[0], first_vline_y[0]]),
            point
        )
        vec2 = utils.getVector(
            np.array([second_vline_x[0], second_vline_y[0]]),
            point
        )

        angle1 = utils.getAngle(vec1, unit_vector)
        angle2 = utils.getAngle(vec2, unit_vector)

        if angle1 <= first_angle and angle2 >= second_angle:
            # If true add index to list marked for deletion
            ind_to_del = np.append(ind_to_del, i)

            if i % 2 == 0:
                # If even add next index
                ind_to_del = np.append(ind_to_del, i+1)

            else:
                # If odd add previous index
                ind_to_del = np.append(ind_to_del, i-1)

    ind_to_del = np.unique(ind_to_del)  # Remove duplicates
    projection_points[ind_to_del] = np.array([0, 0])

    return projection_points


def findProjectionPoints(img, centre_points, nb_points, horn):
    """Find the projection points from the centre point onto the muscle
    layers given the desired number of points.

    Arguments:
    img -- ndarray, image to analyse.
    centre_points -- list[int], coordinates of the centre points (XY).
    nb_points -- int, number of desired projection points, must be a
        multiple of 2.
    horn -- str {left, right}, horn that is being analysed.

    Return:
    projection_points -- ndarray, list of the coordinates of the
        projection points.

    """
    try:
        assert nb_points % 2 == 0

    except AssertionError:
        sys.stderr.write("Error: the number of points needs to be even.\n")
        exit(1)

    angles = np.arange(1, 1 + (nb_points / 2)) * (np.pi / (nb_points / 2))
    projection_points = np.zeros((nb_points * 2, 2), dtype=int)

    if (centre_points[2:4] != np.array([0, 0])).all():
        # The horns are not clearly separated and three points are given
        # Create a vector between the left and right points
        normal = utils.getVector(
            np.array([centre_points[5], centre_points[0]]),
            np.array([centre_points[1], centre_points[4]]),
        )

        # Use the middle point to draw a line
        line_x, _ = separationLine(img.shape, centre_points[2:4], normal)

        for i in range(len(line_x)):
            # Clear half of the image based on the horn
            match horn:
                case "left":
                    img[i, line_x[i]:] = 0

                case "right":
                    img[i, :line_x[i]+1] = 0

    # Select the correct centre point based on the horn
    match horn:
        case "left":
            centre_point = centre_points[0:2]

        case "right":
            centre_point = centre_points[4:6]

    for i, theta in enumerate(angles):
        # Find the line for the given angle
        line_x, line_y = findLineCoordinates(img.shape, centre_point, theta)
        line = img[line_y, line_x]

        # Find the indices of rising and falling edges
        coords = np.where(line[:-1] != line[1:])[0]

        for j in np.arange(0, len(coords), 2):
            coords[j] += 1

        x_coords = line_x[coords]
        y_coords = line_y[coords]

        points = createProjectionPointCoords(x_coords, y_coords, centre_point,
                                             theta)

        projection_points[i * 4: (i + 1) * 4] = points

    if (centre_points[2:4] != np.array([0, 0])).all():
        projection_points = excludeCentralPoints(
            img, centre_points, projection_points, horn)

    return projection_points


def createProjectionPointCoords(x_coords, y_coords, centre_point, theta):
    """Creates the (x, y) pairs of coordinates for the projection points

    The points that are to the right of the centre point are placed first
    in the array. In the case of the vertical line, the top points are
    placed first.

    Arguments:
    x_coords -- list[int], list of coordinates of the projection points on
        the x axis.
    y_coords -- list[int], list of coordinates of the projection points on
        the y axis.
    centre_point -- ndarray, coordinates of the centre point (XY)
    theta -- float, angle at which the projection line is on.

    Return:
    projection_points -- ndarray, list of coordinates of the
        projection points.

    """
    try:
        assert len(x_coords) == len(y_coords)

    except AssertionError:
        sys.stderr.write(
            "Error: x_coords and y_coords should have the same size.\n")
        exit(1)

    transpose = False
    point_list = np.transpose([x_coords, y_coords])

    # Get the indices of points before and after on the first axis
    diff = point_list - centre_point
    neg_indices = np.arange(len(diff))[diff[:, 0] < 0]
    pos_indices = np.arange(len(diff))[diff[:, 0] >= 0]

    if pos_indices.shape[0] <= 1 or neg_indices.shape[0] <= 1:
        # If the x coords are not enough to distinguish use y coords
        point_list = np.flip(point_list, axis=1)  # Flip to [y, x]
        centre_point = np.flip(centre_point)  # Flip to [y, x]

        diff = point_list - centre_point
        neg_indices = np.arange(len(diff))[diff[:, 0] < 0]
        pos_indices = np.arange(len(diff))[diff[:, 0] >= 0]
        transpose = True

    if pos_indices.shape[0] <= 1 or neg_indices.shape[0] <= 1:
        # If some points where not found, assure minimal distance of 1
        point_list = np.flip(point_list, axis=1)  # Flip to [x, y]
        centre_point = np.flip(centre_point)  # Flip to [x, y]

        diff = point_list - centre_point
        neg_indices = np.arange(len(diff))[diff[:, 0] < 0]
        pos_indices = np.arange(len(diff))[diff[:, 0] >= 0]

        if diff[0, 0] < 0:
            # Points need to be added after the centre point
            point_list = np.concatenate(
                (
                    [[centre_point[0] + 2, centre_point[1]]],
                    [[centre_point[0] + 1, centre_point[1]]],
                    point_list,
                )
            )

        elif diff[0, 0] >= 0:
            # Points need to be added before the centre point
            point_list = np.concatenate(
                (
                    point_list,
                    [[centre_point[0] - 1, centre_point[1]]],
                    [[centre_point[0] - 2, centre_point[1]]],
                )
            )

        diff = point_list - centre_point
        neg_indices = np.arange(len(diff))[diff[:, 0] < 0]
        pos_indices = np.arange(len(diff))[diff[:, 0] >= 0]

    # Find the distances of the points from the centre
    distances_neg = np.linalg.norm(diff[neg_indices], axis=1)
    distances_pos = np.linalg.norm(diff[pos_indices], axis=1)

    # Sort indices to ensure that the points are sorted
    # from nearest to furthest
    neg_indices = neg_indices[np.argsort(distances_neg)]
    pos_indices = pos_indices[np.argsort(distances_pos)]

    # Create the sets of points on the inner and outer edges
    if transpose:
        # Flip the points back to [x, y]
        point_list = np.flip(point_list, axis=1)
        centre_point = np.flip(centre_point)

    if theta < np.pi / 2:
        # The positive points are on the left
        first_set = point_list[[neg_indices[0], neg_indices[1]]]
        second_set = point_list[[pos_indices[0], pos_indices[1]]]

    else:
        # The positive points are on the right
        first_set = point_list[[pos_indices[0], pos_indices[1]]]
        second_set = point_list[[neg_indices[0], neg_indices[1]]]

    projection_points = np.concatenate((first_set, second_set))

    return projection_points


def alignBorder(thickness):
    """Organises the thickness array to align the first value
        with the anti-mesometrial border

    Arguments:
    thickness -- ndarray, array of thickness for each angle.

    Return:
    ordered_thickness -- ndarray, ordered thickness array.

    """
    nb_points = len(thickness)

    # Find the two halves
    right_half = np.arange(0, nb_points, 2)
    left_half = np.arange(1, nb_points, 2)

    # Order thickness to go from 0 to 2pi
    ordered_thickness = np.concatenate((thickness[right_half],
                                        thickness[left_half]))

    # Find the maximal thickness which typically is AM border
    max_idx = np.argmax(ordered_thickness)

    if np.isnan(ordered_thickness[max_idx]):
        if ~np.isnan(ordered_thickness[right_half][-1]) and ~np.isnan(
                ordered_thickness[left_half][0]):
            # The nan values are in the middle
            nan_indices = np.where(np.isnan(ordered_thickness))[0]
            max_idx = nan_indices[len(nan_indices) // 2]

        else:
            # The nan values are already on each side of the array
            max_idx = 0

    # Roll array to line up 0 with anti-mesometrial border
    ordered_thickness = np.roll(ordered_thickness, nb_points - max_idx)

    return ordered_thickness


def estimateMuscleThickness(img_stack, centreline, nb_points, slice_nbs, horn):
    """Estimates the muscle thickness of each slice

    Arguments:
    img_stack -- ndarray, stack of images to process (ZXY dimensions).
    centreline -- ndarray, coordinates of the centre point on each slice.
    nb_points -- int, number of points to use for projection.
    slice_nbs -- list[int], number of the slices at which to save
        angular thickness.
    horn -- str {left, right}, horn that is being analysed.

    Return:
    muscle_thickness_array -- ndarray, muscle thickness of each slice.
    slice_thickness_array -- ndarray, muscle thickness at different angles
        for three slices.
    radius_array -- ndarray, average radius of each slice.

    """
    try:
        assert len(img_stack) == len(centreline)

    except AssertionError:
        sys.stderr.write(
            "Error: the image stack and the centreline do not "
            " have the same size.\n"
        )
        exit()

    nb_imgs = len(img_stack)
    muscle_thickness_array = np.zeros(nb_imgs)
    radius_array = np.zeros(nb_imgs)
    slice_thickness_array = list()
    idx_removed_slices = list()

    # Create error log file
    log_filename = "uCT_errors.log"
    logging.basicConfig(
        filename=log_filename,
        filemode="w",
        level=logging.ERROR,
        format="%(levelname)s %(name)s %(message)s",
    )
    logger = logging.getLogger(__name__ + "_" + horn)

    for i, img in enumerate(img_stack):
        try:
            projection_points = findProjectionPoints(
                img, centreline[i, :], nb_points, horn
            )

            diff = np.diff(projection_points, axis=0)
            norm = np.linalg.norm(diff, axis=1)
            thickness = norm[np.arange(0, projection_points.shape[0], 2)]

            exclusion_indices = np.where(thickness == 0)[0]
            thickness[exclusion_indices] = float('nan')
            muscle_thickness_array[i] = np.mean(
                thickness[~np.isnan(thickness)])

            # Get points in the inner layer of muscles
            inner_points = projection_points[
                np.arange(0, projection_points.shape[0], 2)
            ]
            # Remove projection points in the exclusion zone
            inner_indices = np.where(inner_points != np.array([0, 0]))[0]
            inner_points = inner_points[inner_indices]

            # Select the correct centre point based on the horn
            match horn:
                case "left":
                    centre_point = centreline[i, 0:2]

                case "right":
                    centre_point = centreline[i, 4:6]

            # Estimate average radius for the slice
            radius_array[i] = np.mean(
                np.linalg.norm(inner_points - centre_point, axis=1)
            )

            if i in slice_nbs:
                ordered_thickness = alignBorder(thickness)
                slice_thickness_array.append(ordered_thickness)

        except Exception as err:
            # Write to stderr and add slice number to the remove list
            sys.stderr.write(
                "Warning: unable to process image number {}\n".format(i))
            idx_removed_slices.append(i)
            logger.exception(err)
    # Remove the slices the values of the slices that were not processed
    muscle_thickness_array = np.delete(muscle_thickness_array,
                                       idx_removed_slices)

    radius_array = np.delete(radius_array, idx_removed_slices)
    return muscle_thickness_array, np.transpose(
        slice_thickness_array), radius_array
