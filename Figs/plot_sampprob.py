import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib
from astropy.table import Table

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # data
    gor22 = Table.read("data/edf_fits_param_edit.csv")
    gor22["b3"] = gor22["c3"] / np.square(gor22["gamma"])
    gor22["av"] = gor22["ebv"] * gor22["rv"]
    gor22["av_func"] = gor22["av"] * np.sqrt(
        np.square(gor22["ebv_func"] / gor22["ebv"])
        + np.square(gor22["rv_func"] / gor22["rv"])
    )

    # plots
    fontsize = 14

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    figsize = (5.5, 5.0)
    fig, ax = pyplot.subplots(figsize=figsize)

    # R(V) versus A(V)
    ax.errorbar(
        gor22["rv"], gor22["av"], fmt="go", xerr=gor22["rv_func"], yerr=gor22["av_func"]
    )
    ax.set_xlabel(r"$R(V)$")
    ax.set_ylabel(r"$A(V)$")
    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")

    # ax.legend()

    fig.tight_layout()

    outbase = "_sampprop"
    if args.png:
        fig.savefig(outbase + ".png")
    elif args.pdf:
        fig.savefig(outbase + ".pdf")
    else:
        pyplot.show()
