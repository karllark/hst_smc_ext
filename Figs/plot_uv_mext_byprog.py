#!/usr/bin/env python

# import pkg_resources
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec

# import astropy.units as u

from plot_extinction import plot_ext_stack

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--forecor", help="use the foreground corrected curves", action="store_true")
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

    datapath = "/home/kgordon/Python/hst_smc_ext/fits/"

    plot_ext_stack(
        "data/smc_ext_iue.dat",
        ax[0],
        locpath=datapath,
        fontsize=fontsize,
        forecor=args.forecor,
    )
    ax[0].text(
        0.7, 0.05, "IUE", fontsize=1.5 * fontsize, transform=ax[0].transAxes, bbox=props
    )
    ax[0].set_ylim(-5.0, 65.0)

    plot_ext_stack(
        "data/smc_ext_stis1.dat",
        ax[1],
        locpath=datapath,
        fontsize=fontsize,
        forecor=args.forecor,
    )
    ax[1].text(
        0.7, 0.06, "STIS1", fontsize=1.5 * fontsize, transform=ax[1].transAxes, bbox=props
    )
    ax[1].set_ylim(-5.0, 45.0)

    plot_ext_stack(
        "data/smc_ext_stis2.dat",
        ax[2],
        locpath=datapath,
        fontsize=fontsize,
        forecor=args.forecor,
    )
    ax[2].text(
        0.7, 0.05, "STIS2", fontsize=1.5 * fontsize, transform=ax[2].transAxes, bbox=props
    )
    ax[2].set_ylim(-5.0, 80.0)

    plot_ext_stack(
        "data/smc_ext_stis3.dat",
        ax[3],
        locpath=datapath,
        fontsize=fontsize,
        forecor=args.forecor,
    )
    ax[3].text(
        0.7, 0.08, "STIS3", fontsize=1.5 * fontsize, transform=ax[3].transAxes, bbox=props
    )
    ax[3].set_ylim(-5.0, 25.0)

    for i in range(4):
        ax[i].set_xlim(1. / 0.7, 1. / 0.11)
    #     ylim = ax[i].get_ylim()
        # ax[i].set_xscale("log")
        # ax[i].set_yscale("log")
        # ax[i].set_xlabel(r"$\lambda$ [$\mu m$]", fontsize=1.3 * fontsize)
    #     # )
        ax[i].tick_params("both", length=10, width=2, which="major")
        ax[i].tick_params("both", length=5, width=1, which="minor")
    #
    ax[0].set_ylabel(
        r"$E(\lambda - V)/E(B-V)$ + offset",
        fontsize=1.3 * fontsize,
    )
    #
    # ax[2].set_ylabel(
    #     r"$F(\lambda)/F(2600~\AA)$ + offset [$ergs\ cm^{-2}\ s\ \AA$]",
    #     fontsize=1.3 * fontsize,
    # )
    #
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
    save_str = "smc_stars_ext_byprog"
    if args.forecor:
        save_str = f"{save_str}_forecor"
    if args.png:
        fig.savefig(f"{save_str}.png")
    elif args.pdf:
        fig.savefig(f"{save_str}.pdf")
    else:
        plt.show()
