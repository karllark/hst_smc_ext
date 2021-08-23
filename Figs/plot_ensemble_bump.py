import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from astropy.table import Table

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    locpath = ""

    fontsize = 10

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    # get the new results
    gor22 = Table.read("data/edf_fits_param_edit.csv")
    gor22["b3"] = gor22["c3"] / np.square(gor22["gamma"])
    gor22["av"] = gor22["ebv"] * gor22["rv"]

    gor22["abump"] = gor22["b3"] / gor22["rv"]

    gor22["nhiav"] = (10 ** gor22["logN(HI)"]) / gor22["av"]

    # get the IR information
    irtab = Table.read(
        "data/forgordon_fit_param_smc.txt", format="ascii.commented_header"
    )
    # get translation of names
    ntrans = Table.read(
        "data/smc_stars_all_ed_to_extstar.dat", format="ascii.commented_header"
    )

    tpah = np.full((len(gor22)), -1.0)
    for i, cname in enumerate(ntrans["name2"]):
        if cname.upper() in irtab["Star"]:
            mindx = np.where(ntrans["name1"][i].upper() == gor22["name"])
            nindx = np.where(cname.upper() == irtab["Star"])
            tpah[mindx] = irtab["q_PAH"][nindx]

    gor22["qpah"] = tpah

    fig, ax = plt.subplots(nrows=2, ncols=2, sharex=False, figsize=(10, 10))

    # plot bump versus nhiav
    g22size = 10
    px = 0
    py = 0
    ax[px, py].plot(
        gor22["abump"],
        gor22["nhiav"],
        "*",
        color="tab:green",
        label="this work SMC",
        alpha=0.75,
        markersize=g22size,
    )
    ax[px, py].set_xlabel("A(bump)/AV", fontsize=1.3 * fontsize)
    ax[px, py].set_ylabel("N(HI)/AV", fontsize=1.3 * fontsize)

    # plot bump versus c2/rv
    g22size = 10
    px = 0
    py = 1
    ax[px, py].plot(
        gor22["abump"],
        gor22["c2"] / gor22["rv"],
        "*",
        color="tab:green",
        label="this work SMC",
        alpha=0.75,
        markersize=g22size,
    )
    ax[px, py].set_xlabel("A(bump)/AV", fontsize=1.3 * fontsize)
    ax[px, py].set_ylabel("C2/RV", fontsize=1.3 * fontsize)

    # plot bump versus qpah
    g22size = 10
    px = 1
    py = 0
    gvals = gor22["qpah"] > -0.9
    ax[px, py].plot(
        gor22["abump"][gvals],
        gor22["qpah"][gvals],
        "*",
        color="tab:green",
        label="this work SMC",
        alpha=0.75,
        markersize=g22size,
    )
    ax[px, py].set_xlabel("A(bump)/AV", fontsize=1.3 * fontsize)
    ax[px, py].set_ylabel(r"$q_{PAH}$", fontsize=1.3 * fontsize)

    # plot bump versus b3
    g22size = 10
    px = 1
    py = 1
    gvals = gor22["qpah"] > -0.9
    ax[px, py].plot(
        gor22["b3"][gvals],
        gor22["qpah"][gvals],
        "*",
        color="tab:green",
        label="this work SMC",
        alpha=0.75,
        markersize=g22size,
    )
    ax[px, py].set_xlabel("B3", fontsize=1.3 * fontsize)
    ax[px, py].set_ylabel(r"$q_{PAH}$", fontsize=1.3 * fontsize)

    plt.tight_layout()

    outbase = "ensemble_bump"
    if args.png:
        fig.savefig(outbase + ".png")
    elif args.pdf:
        fig.savefig(outbase + ".pdf")
    else:
        plt.show()
