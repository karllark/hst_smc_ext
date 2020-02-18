import numpy as np
import matplotlib.pyplot as pyplot
import argparse
import matplotlib

from extdata import *


def plot_ext_stack(filelist, ax, idlsave=False, locpath="./", fontsize=14):

    f = open(locpath + filelist, "r")
    file_lines = list(f)
    starnames = []
    extdata = []
    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            starnames.append(name)
            text = ExtData()
            if idlsave:
                text.read_ext_data_idlsave(
                    locpath + starnames[-1].upper() + "_EXTCURVE.save"
#                    locpath + "FITSPEC_" + starnames[-1].upper() + "_ATLAS.save"
                )
            else:
                text.read_ext_data(locpath + starnames[-1] + "_ext.fits")
            extdata.append(text)

    norm_wave_range = [0.4, 0.6]
    ann_wave_range = 1.0 / np.array([8.65, 8.45])
    col_vals = ["b", "g", "r", "m", "c", "y"]
    offset_val = 5.0

    n_stars = len(starnames)
    for i in range(n_stars):
        k = i

        if idlsave:
            norm_val = 1.0
        else:
            norm_indxs = np.where(
                (extdata[k].ext_waves["BANDS"] >= norm_wave_range[0])
                & (extdata[k].ext_waves["BANDS"] <= norm_wave_range[1])
            )
            norm_val = 1.0 / np.average(extdata[k].ext_curve["BANDS"][norm_indxs])

            ax.plot(
                1.0 / extdata[k].ext_waves["BANDS"],
                extdata[k].ext_curve["BANDS"] * norm_val + offset_val * i,
                col_vals[i % 6] + "o",
            )

        ax.plot(
            1.0 / extdata[k].ext_waves["STIS"],
            extdata[k].ext_curve["STIS"] * norm_val + offset_val * i,
            col_vals[i % 6] + "-",
        )

        if idlsave:
            indxs, = np.where(1.0 / extdata[k].ext_waves["MODEL"] < 9.0)
            ax.plot(
                1.0 / extdata[k].ext_waves["MODEL"][indxs],
                extdata[k].ext_curve["MODEL"][indxs] * norm_val + offset_val * i,
                col_vals[i % 6] + "--",
            )

        ann_indxs = np.where(
            (extdata[k].ext_waves["STIS"] >= ann_wave_range[0])
            & (extdata[k].ext_waves["STIS"] <= ann_wave_range[1])
        )
        ann_val = np.average(extdata[k].ext_curve["STIS"][ann_indxs])
        ann_val *= norm_val
        ann_val += offset_val * i
        ax.annotate(
            starnames[k],
            xy=(8.8, ann_val),
            xytext=(9.3, ann_val),
            verticalalignment="center",
            arrowprops=dict(facecolor=col_vals[i % 6], shrink=0.1),
            fontsize=0.85 * fontsize,
            rotation=-45.0,
        )

    ax.set_xlim(1.0 / 2.5, 10.5)
    ax.set_ylim(-5.0, 20.0 + offset_val * n_stars)


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of stars to plot")
    parser.add_argument(
        "-i",
        "--idlsave",
        help="extinction data in an IDL save file",
        action="store_true",
    )
    parser.add_argument(
        "-p", "--png", help="save figure as a png file", action="store_true"
    )
    parser.add_argument(
        "-e", "--eps", help="save figure as an eps file", action="store_true"
    )
    parser.add_argument("--path", help="path to filelist and extinction curves")
    args = parser.parse_args()

    if args.path:
        locpath = args.path + "/"
    else:
        locpath = ""

    fontsize = 18

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    fig, ax = pyplot.subplots(nrows=1, ncols=1, figsize=(10, 12))

    plot_ext_stack(args.filelist, ax, idlsave=args.idlsave, locpath=locpath)

    ax.set_xlabel("$1/\lambda$ [$\mu m^{-1}$]", fontsize=1.3 * fontsize)
    ax.set_ylabel("$E(\lambda - V)/E(B - V)$ + offset", fontsize=1.3 * fontsize)

    pyplot.show()
