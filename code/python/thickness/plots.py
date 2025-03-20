#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plots.py

Plot functions for the thickness package
Author: Mathias Roesler
Date: 03/25
"""

import matplotlib.pyplot as plt
import numpy as np

from thickness.constants import LEFT, BOTTOM, RIGHT, COLOURS, Y_LABELS


def plotProjectionPoints(img, centre, projection_points):
    """Plots the centre point and projection points of an image

    Arguments:
    img -- ndarray, image to display.
    centre -- ndarray, coordinates of the centre point (XY).
    projection_points -- ndarray, coordinates of the projection points (XY).

    Return:

    """
    fig, ax = plt.subplots(dpi=300)

    plt.imshow(img, cmap="gray")  # Plot image
    ax.plot(centre[0], centre[1], ".r")  # Plot centre

    for point in projection_points:
        ax.plot(point[0], point[1], ".b")

    plt.show()


def plotMuscleThickness(muscle_thickness, errors):
    """Plots the muscle thickness of both horns on the same plot

    The number of points is normalised so that both sets are shown
    between 0 and 1

    Arguments:
    muscle_thickness -- dict(ndarray), thickness of the horns.
    errors -- dict(ndarray), errors for the muscle thickness.

    Return:

    """
    fig, ax = plt.subplots(dpi=300)
    colors = {"left": "tab:blue", "right": "tab:red"}

    for horn in muscle_thickness.keys():
        horn_thickness = muscle_thickness[horn]
        error_bars = errors[horn]
        horn_length = len(horn_thickness)

        ax.errorbar(
            np.linspace(0, 1, horn_length),
            horn_thickness,
            yerr=error_bars,
            linewidth=2,
            label="{} horn".format(horn.capitalize()),
            color=colors[horn],
        )

    # Reset x-axis ticks
    plt.xticks(
        ticks=[0, 0.25, 0.6, 1],
        labels=["Cervix", "Cervical end", "Centre", "Ovarian end"],
    )

    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.xlabel("Locations")
    plt.ylabel("Muscle thickness (mm)")
    plt.title("Average thickness in the uterus")
    plt.legend()

    plt.subplots_adjust(left=LEFT, right=RIGHT, bottom=BOTTOM)
    plt.show()


def plotAngularThickness(slice_thickness, projection=False, uCT_flag=True):
    """Plots the muscle thickness of one slice as a function of the
    angle theta

    Arguments:
    slice_thickness -- dict(ndarray), array containing the angluar thickness
        of four slices for each horn.
    projection -- str, projection type of the plot, default value False.
    uCT_flag -- bool, True if the data is from microCT False if the data is
        from histology, default value True.

    Return:

    """
    fig, ax = plt.subplots(
        len(slice_thickness.keys()),
        1,
        subplot_kw={"polar": projection},
        dpi=300,
    )

    if not hasattr(ax, "__len__"):
        # If only one subplot is created
        ax = [ax]  # Convert to list for rest of code to work

    if uCT_flag:
        data_type = r"$\mu$CT"

    else:
        data_type = "histology"

    for i, horn in enumerate(slice_thickness.keys()):
        y_values = slice_thickness[horn]

        # Create x-axis values so that everything is normalised
        nb_points = y_values.shape[0]
        x_values = np.linspace(0, 2 * np.pi, nb_points, endpoint=False)

        ax[i].plot(
            x_values,
            y_values[:, 0],
            linestyle="dashdot",
            color="tab:gray",
            label="Cervix",
            linewidth=2,
        )
        ax[i].plot(
            x_values,
            y_values[:, 1],
            linestyle="solid",
            color="tab:orange",
            label="Cervical end",
            linewidth=2,
        )
        ax[i].plot(
            x_values,
            y_values[:, 2],
            linestyle="dashed",
            color="tab:purple",
            label="Centre",
            linewidth=2,
        )
        ax[i].plot(
            x_values,
            y_values[:, 3],
            linestyle="dotted",
            color="tab:green",
            label="Ovarian end",
            linewidth=2,
        )
        ax[i].set_title(
            r"{} horn thickness (from {} data)".format(
                horn.capitalize(),
                data_type,
            )
        )

        if projection:
            ax[i].set_rlabel_position(-22.5)  # Move radial labels

            ax[i].set_rmax(1.1)  # Set radial max
            ticks = plt.xticks()[0]

            # Set labels and legends
            angle = np.deg2rad(25)
            plt.legend(
                loc="lower left",
                bbox_to_anchor=(0.5 + np.cos(angle) / 2, 0.5 + np.sin(angle) / 2),
            )

            plt.xticks(
                ticks=ticks,
                labels=[
                    "0",
                    r"$\frac{\pi}{4}$",
                    r"$\frac{\pi}{2}$",
                    r"$\frac{3\pi}{4}$",
                    r"$\pi$",
                    r"$\frac{5\pi}{4}$",
                    r"$\frac{3\pi}{2}$",
                    r"$\frac{7\pi}{4}$",
                ],
            )

        else:
            plt.xlim([0, 2 * np.pi])
            plt.ylim([0, 1.1])
            ticks = np.linspace(0, 2 * np.pi, 9)
            plt.legend(loc="upper left")

            plt.xticks(
                ticks=ticks,
                labels=[
                    "0",
                    r"$\frac{\pi}{4}$",
                    r"$\frac{\pi}{2}$",
                    r"$\frac{3\pi}{4}$",
                    r"$\pi$",
                    r"$\frac{5\pi}{4}$",
                    r"$\frac{3\pi}{2}$",
                    r"$\frac{7\pi}{4}$",
                    r"2$\pi$",
                ],
            )

            plt.ylabel("Muscle thickness (mm)")
            plt.xlabel(r"Angle $\theta$ (rad)")

    plt.subplots_adjust(left=LEFT, right=RIGHT, bottom=BOTTOM)
    plt.show()


def plotData(data, metric):
    """Plots the selected data.

    Arguments:
    data -- dict(list(float))), dictionnary with estrus phases as keys and
    lists of metric values as values

    Return:

    """
    fig, ax = plt.subplots(dpi=300)

    for i, stage in enumerate(data.keys()):
        nb_samples = len(data[stage])
        np.random.seed(12)  # Reset random seed for all stages to be identical
        jitter = np.random.uniform(-0.1, 0.1, nb_samples)

        plt.errorbar(
            (i + 1) * np.ones(nb_samples) + jitter,
            data[stage][:, 0],
            data[stage][:, 1],
            c=COLOURS[stage],
            marker=".",
            linestyle="",
            capsize=3.0,
        )

    # Reset x-axis ticks
    plt.xticks(
        ticks=[1, 2, 3, 4],
        labels=[estrus.capitalize() for estrus in data.keys()],
    )
    plt.xlim([0.5, 4.5])
    plt.ylim(bottom=0)

    plt.ylabel(Y_LABELS[metric])
    plt.subplots_adjust(left=LEFT, right=RIGHT, bottom=BOTTOM)
    plt.show()
