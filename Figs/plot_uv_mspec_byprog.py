#!/usr/bin/env python

# import pkg_resources
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec

# import numpy as np

import astropy.units as u

from measure_extinction.stardata import StarData
from measure_extinction.utils.helpers import get_full_starfile

from plot_uv_mspec import get_starnames


def plot_sspec(starnames, ax):
    # plot the spectra in two columns
    col_vals = ["b", "g", "c"]
    n_cols = len(col_vals)
    for k, cstarname in enumerate(starnames):
        fstarname, file_path = get_full_starfile(cstarname)
        starobs = StarData(fstarname, path=file_path)

        if "IUE" not in starobs.data.keys():
            bkey = "STIS"
        else:
            bkey = "IUE"

        yoff = 2.5**k

        # remove negative values as they look horrible on stacked log plots
        nvals = starobs.data[bkey].fluxes < 0
        starobs.data[bkey].npts[nvals] = 0

        if cstarname in ["mr12-star09", "mr12-star10", "mr12-star11"]:
            starobs.data[bkey].rebin_constres([0.1, 0.34] * u.micron, 500.0)

        starobs.plot(
            ax,
            norm_wave_range=[0.2, 0.3] * u.micron,
            yoffset=yoff,
            pcolor=col_vals[k % n_cols],
            annotate_key=bkey,
            annotate_wave_range=[0.25, 0.27] * u.micron,
            annotate_text=cstarname,
            annotate_rotation=-10.0,
            annotate_yoffset=0.0,
            fontsize=12,
        )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--path", help="path to star files", default="DAT_Files/")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # plotting setup for easier to read plots
    fontsize = 14
    font = {"size": fontsize}
    mpl.rc("font", **font)
    mpl.rc("lines", linewidth=1.5)
    mpl.rc("axes", linewidth=3)
    mpl.rc("xtick.major", width=3)
    mpl.rc("xtick.minor", width=3)
    mpl.rc("ytick.major", width=3)
    mpl.rc("ytick.minor", width=3)

    fig, ax = plt.subplots(figsize=(18, 12))
    gs = gridspec.GridSpec(20, 2)
    ax = []
    ax.append(plt.subplot(gs[0:12, 0]))
    ax.append(plt.subplot(gs[12:20, 0]))
    ax.append(plt.subplot(gs[0:15, 1]))
    ax.append(plt.subplot(gs[15:20, 1]))

    # fig, max = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    # ax = [max[0, 0], max[0, 1], max[1, 0], max[0, 0]]

    props = dict(boxstyle="round", facecolor="lightgrey", alpha=0.5)

    starnames = get_starnames("data/smc_stars_iue.dat")
    plot_sspec(starnames, ax[0])
    ax[0].text(
        0.7, 0.05, "IUE", fontsize=1.5 * fontsize, transform=ax[0].transAxes, bbox=props
    )

    starnames = get_starnames("data/smc_stars_stis1.dat")
    plot_sspec(starnames, ax[1])
    ax[1].text(
        0.7, 0.06, "STIS1", fontsize=1.5 * fontsize, transform=ax[1].transAxes, bbox=props
    )

    starnames = get_starnames("data/smc_stars_stis2.dat")
    plot_sspec(starnames, ax[2])
    ax[2].text(
        0.7, 0.05, "STIS2", fontsize=1.5 * fontsize, transform=ax[2].transAxes, bbox=props
    )

    starnames = get_starnames("data/smc_stars_stis3.dat")
    plot_sspec(starnames, ax[3])
    ax[3].text(
        0.7, 0.08, "STIS3", fontsize=1.5 * fontsize, transform=ax[3].transAxes, bbox=props
    )

    for i in range(4):
        ax[i].set_xlim(0.1, 0.325)
        ylim = ax[i].get_ylim()
        ax[i].set_ylim(1e-1, ylim[1])
        ax[i].set_xscale("linear")
        ax[i].set_yscale("log")
        ax[i].set_xlabel(r"$\lambda$ [$\mu m$]", fontsize=1.3 * fontsize)
        # )
        ax[i].tick_params("both", length=10, width=2, which="major")
        ax[i].tick_params("both", length=5, width=1, which="minor")

    ax[0].set_ylabel(
        r"$F(\lambda)/F(2600~\AA)$ + offset [$ergs\ cm^{-2}\ s\ \AA$]",
        fontsize=1.3 * fontsize,
    )

    ax[2].set_ylabel(
        r"$F(\lambda)/F(2600~\AA)$ + offset [$ergs\ cm^{-2}\ s\ \AA$]",
        fontsize=1.3 * fontsize,
    )

    ax[0].xaxis.tick_top()
    ax[0].xaxis.set_label_position("top")

    ax[2].xaxis.tick_top()
    ax[2].xaxis.set_label_position("top")

    ax[2].yaxis.tick_right()
    ax[2].yaxis.set_label_position("right")

    ax[3].yaxis.tick_right()
    ax[3].yaxis.set_label_position("right")

    # use the whitespace better
    fig.tight_layout()

    # plot or save to a file
    save_str = "smc_stars_spec_byprog"
    if args.png:
        fig.savefig(f"{save_str}.png")
    elif args.pdf:
        fig.savefig(f"{save_str}.pdf")
    else:
        plt.show()
