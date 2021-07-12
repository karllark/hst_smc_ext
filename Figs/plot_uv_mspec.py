#!/usr/bin/env python

from __future__ import absolute_import, division, print_function, unicode_literals

# import pkg_resources
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

import astropy.units as u

from measure_extinction.stardata import StarData
from measure_extinction.utils.helpers import get_full_starfile


def plot_mspec_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of stars to plot")
    parser.add_argument("--path", help="path to star files", default="./")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    return parser


if __name__ == "__main__":

    # commandline parser
    parser = plot_mspec_parser()
    args = parser.parse_args()

    # get the names of the stars
    f = open(args.filelist, "r")
    file_lines = list(f)
    starnames = []
    spslopes = []
    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            starnames.append(name)

            # use the IUE spectra to make a measure of the UV spectral slope
            # used to sort for a nicer plot
            fstarname, file_path = get_full_starfile(name)
            starobs = StarData(fstarname, path=file_path)

            gvals1 = np.logical_and(
                starobs.data["IUE"].waves > 0.15 * u.micron,
                starobs.data["IUE"].waves < 0.20 * u.micron,
            )
            gavls1 = np.logical_and(gvals1, starobs.data["IUE"].npts > 0)
            gvals2 = np.logical_and(
                starobs.data["IUE"].waves > 0.25 * u.micron,
                starobs.data["IUE"].waves < 0.30 * u.micron,
            )
            gavls2 = np.logical_and(gvals2, starobs.data["IUE"].npts > 0)
            spslopes.append(
                np.median(starobs.data["IUE"].fluxes[gvals1])
                / np.median(starobs.data["IUE"].fluxes[gvals2])
            )

    # sort starnames by UV slope
    starnames = np.array(starnames)
    starnames = starnames[np.argsort(spslopes)]

    # plotting setup for easier to read plots
    fontsize = 14
    font = {"size": fontsize}
    mpl.rc("font", **font)
    mpl.rc("lines", linewidth=1)
    mpl.rc("axes", linewidth=2)
    mpl.rc("xtick.major", width=2)
    mpl.rc("xtick.minor", width=2)
    mpl.rc("ytick.major", width=2)
    mpl.rc("ytick.minor", width=2)

    # setup the plot
    fig, ax = plt.subplots(ncols=2, figsize=(13, 10))

    # plot the spectra in two columns
    half_num = len(starnames) // 2 + 1
    col_vals = ["b", "g", "c"]
    n_cols = len(col_vals)
    for k, cstarname in enumerate(starnames):
        fstarname, file_path = get_full_starfile(cstarname)
        starobs = StarData(fstarname, path=file_path)

        # remove negative values as they look horrible on stacked log plots
        nvals = starobs.data["IUE"].fluxes < 0
        starobs.data["IUE"].npts[nvals] = 0

        if k // half_num > 0:
            yoff = 2.5 ** (k - half_num)
        else:
            yoff = 2.5 ** k
        starobs.plot(
            ax[k // half_num],
            norm_wave_range=[0.2, 0.3] * u.micron,
            yoffset=yoff,
            pcolor=col_vals[k % n_cols],
            annotate_key="IUE",
            annotate_wave_range=[0.25, 0.27] * u.micron,
            annotate_text=cstarname,
            annotate_rotation=-10.0,
            annotate_yoffset=0.0,
            fontsize=12,
        )

    # finish configuring the plot
    for i in range(2):
        ax[i].set_xlim(0.1, 0.325)
        ylim = ax[i].get_ylim()
        ax[i].set_ylim(1e-1, ylim[1])
        ax[i].set_xscale("linear")
        ax[i].set_yscale("log")
        ax[i].set_xlabel(r"$\lambda$ [$\mu m$]", fontsize=1.3 * fontsize)
        ax[i].set_ylabel(
            r"$F(\lambda)/F(2500~\AA)$ + offset [$ergs\ cm^{-2}\ s\ \AA$]",
            fontsize=1.3 * fontsize,
        )
        ax[i].tick_params("both", length=10, width=2, which="major")
        ax[i].tick_params("both", length=5, width=1, which="minor")

    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")

    # use the whitespace better
    fig.tight_layout()

    # plot or save to a file
    save_str = "Figs/smcext_uvspec"
    if args.png:
        fig.savefig(f"{save_str}.png")
    elif args.pdf:
        fig.savefig(f"{save_str}.pdf")
    else:
        plt.show()
