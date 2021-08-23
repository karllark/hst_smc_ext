# Program to plot the ensemble properties of dust extinction fits
#  - FM parameters

import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from astropy.table import Table

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--av", help="plot versus A(V)", action="store_true")
    parser.add_argument("--rv", help="plot versus R(V)", action="store_true")
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
    print(gor22.colnames)
    gor22["b3"] = gor22["c3"] / np.square(gor22["gamma"])
    gor22["av"] = gor22["ebv"] * gor22["rv"]

    # get the Gordon et al. 2009 results
    gor09_dust = Table.read(
        "litdata/Gordon09/gor09_table2.dat", format="ascii.commented_header"
    )
    gor09_fm = Table.read(
        "litdata/Gordon09/gor09_table3.dat", format="ascii.commented_header"
    )

    gor09 = Table()
    gor09["rv"] = np.array(gor09_dust["rv"])
    gor09["av"] = np.array(gor09_dust["av"])
    gor09["ebv"] = np.array(gor09_dust["ebv"])
    gor09["c1"] = (np.array(gor09_fm["cav1"]) - 1.0) * gor09["rv"]
    gor09["c2"] = np.array(gor09_fm["cav2"]) * gor09["rv"]
    gor09["c3"] = np.array(gor09_fm["cav3"]) * gor09["rv"]
    gor09["c4"] = np.array(gor09_fm["cav4"]) * gor09["rv"]
    gor09["x0"] = np.array(gor09_fm["x0"])
    gor09["gamma"] = np.array(gor09_fm["gamma"])
    gor09["b3"] = gor09["c3"] / np.square(gor09["gamma"])

    # get Fitzpatrick & Massa 2007 results
    fm07_dust = Table.read("litdata/Fitzpatrick07/fm07_table3.dat", format="ascii.cds")
    fm07_fm = Table.read("litdata/Fitzpatrick07/fm07_table4.dat", format="ascii.cds")

    fm07 = Table()
    fm07["rv"] = fm07_fm["R(V)"]
    fm07["ebv"] = fm07_dust["E(B-V)"]
    fm07["av"] = fm07["ebv"] * fm07["rv"]
    fm07["c1"] = fm07_fm["c1"]
    fm07["c2"] = fm07_fm["c2"]
    fm07["c3"] = fm07_fm["c3"]
    fm07["c4"] = fm07_fm["c4"]
    fm07["c5"] = fm07_fm["c5"]
    fm07["x0"] = fm07_fm["x0"]
    fm07["gamma"] = fm07_fm["gamma"]
    fm07["b3"] = fm07["c3"] / np.square(fm07["gamma"])

    # get the Valencic04 results
    val04_dust = Table.read(
        "litdata/Valencic04/valencic04_datafile4.txt", format="ascii.cds"
    )
    val04_fm = Table.read(
        "litdata/Valencic04/valencic04_datafile5.txt", format="ascii.cds"
    )

    val04 = Table()
    val04["rv"] = np.array(val04_dust["R(V)"])
    val04["av"] = np.array(val04_dust["A(V)"])
    val04["ebv"] = np.array(val04_dust["E(B-V)"])
    val04["c1"] = (np.array(val04_fm["c1/R(V)+1.0"]) - 1.0) * val04["rv"]
    val04["c2"] = np.array(val04_fm["c2/R(V)"]) * val04["rv"]
    val04["c3"] = np.array(val04_fm["c3/R(V)"]) * val04["rv"]
    val04["c4"] = np.array(val04_fm["c4/R(V)"]) * val04["rv"]
    val04["x0"] = np.array(val04_fm["x0"])
    val04["gamma"] = np.array(val04_fm["gamma"])
    val04["b3"] = val04["c3"] / np.square(val04["gamma"])

    # get the Gordon03 results
    gor03_smc = Table.read(
        "litdata/Gordon03/gordon03_smc.dat",
        format="ascii.commented_header",
        header_start=0,
    )
    gor03_smc["b3"] = gor03_smc["c3"] / np.square(gor03_smc["gamma"])
    gor03_smc["av"] = gor03_smc["ebv"] * gor03_smc["rv"]
    gor03_lmc = Table.read(
        "litdata/Gordon03/gordon03_lmc.dat",
        format="ascii.commented_header",
        header_start=0,
    )
    gor03_lmc["b3"] = gor03_lmc["c3"] / np.square(gor03_lmc["gamma"])
    gor03_lmc["av"] = gor03_lmc["ebv"] * gor03_lmc["rv"]

    fig, ax = plt.subplots(nrows=2, ncols=4, sharex=False, figsize=(16, 8))

    mwsize = 4
    g03size = 6
    g22size = 10
    ploty_labels = [
        "$C_1$",
        "$C_2$",
        "$C_3$",
        r"$B_3 = C_3 / \gamma^2$",
        "$C_4$",
        "$C_5$",
        "$x_o$",
        r"$\gamma$",
    ]
    ploty_tags = ["c1", "c2", "c3", "b3", "c4", "c5", "x0", "gamma"]
    pi = [0, 1, 2, 3, 4, 5, 6, 7]
    for i in range(8):
        px, py = divmod(pi[i], 4)
        if args.rv:
            pxcol = "rv"
            pxlabel = "$R(V)$"
        elif args.av:
            pxcol = "av"
            pxlabel = "$A(V)$"
        else:
            pxcol = "ebv"
            pxlabel = "$E(B-V)$"

        if ploty_tags[i] != "c5":
            ax[px, py].plot(
                val04[pxcol],
                val04[ploty_tags[i]],
                "^",
                color="tab:blue",
                label="V04 MW",
                alpha=0.25,
                markersize=mwsize,
            )
            ax[px, py].plot(
                gor09[pxcol],
                gor09[ploty_tags[i]],
                "v",
                color="tab:orange",
                label="GCC09 MW",
                alpha=0.25,
                markersize=mwsize,
            )
            ax[px, py].plot(
                gor03_lmc[pxcol],
                gor03_lmc[ploty_tags[i]],
                "s",
                color="tab:cyan",
                label="G03 LMC",
                alpha=0.5,
                markersize=g03size,
            )
            ax[px, py].plot(
                gor03_smc[pxcol],
                gor03_smc[ploty_tags[i]],
                "o",
                color="tab:pink",
                label="G03 SMC",
                alpha=0.5,
                markersize=g03size,
            )

        ax[px, py].plot(
            fm07[pxcol],
            fm07[ploty_tags[i]],
            ">",
            label="FM07 MW",
            alpha=0.25,
            markersize=mwsize,
            color="tab:purple",
        )
        ax[px, py].plot(
            gor22[pxcol],
            gor22[ploty_tags[i]],
            "*",
            color="tab:green",
            label="this work SMC",
            alpha=0.75,
            markersize=g22size,
        )

        ax[px, py].set_xlabel(pxlabel, fontsize=1.3 * fontsize)
        ax[px, py].set_ylabel(ploty_labels[i], fontsize=1.3 * fontsize)

        if i == 1:
            ax[px, py].legend()

    if args.rv:
        save_str = "_rv"
    else:
        save_str = "_ebv"

    plt.tight_layout()

    outbase = "ensemble_props_fm90"
    if args.rv:
        outbase += "_rv"
    elif args.av:
        outbase += "_av"
    else:
        outbase += "_ebv"
    if args.png:
        fig.savefig(outbase + ".png")
    elif args.pdf:
        fig.savefig(outbase + ".pdf")
    else:
        plt.show()
