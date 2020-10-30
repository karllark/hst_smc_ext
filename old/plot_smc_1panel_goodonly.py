import matplotlib.pyplot as pyplot
import argparse
import matplotlib
import numpy as np
from extdata import ExtData
from matplotlib.ticker import ScalarFormatter


def plot_ext_stack(filelist, ax, idlsave=False, locpath="./", fontsize=14):

    f = open(filelist, "r")
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
                    locpath
                    + starnames[-1].upper()
                    + "_EXTCURVE.save"
                    #                    locpath + "FITSPEC_" + starnames[-1].upper() + "_ATLAS.save"
                )
            else:
                text.read_ext_data(locpath + starnames[-1] + "_ext.fits")
            extdata.append(text)

    norm_wave_range = [0.4, 0.6]
    col_vals = ["b", "g", "r", "m", "c", "y"]
    ncol = len(col_vals)
    offset_val = 5.0

    n_stars = len(starnames)
    for i in range(n_stars):
        k = i

        norm_indxs = np.where(
            (extdata[k].ext_waves["BANDS"] >= norm_wave_range[0])
            & (extdata[k].ext_waves["BANDS"] <= norm_wave_range[1])
        )
        norm_val = 1.0 / np.average(extdata[k].ext_curve["BANDS"][norm_indxs])

        norm_val = 1.0
        # print(norm_val)

        ax.plot(
            extdata[k].ext_waves["BANDS"] * 1e3,
            extdata[k].ext_curve["BANDS"] * norm_val + offset_val * i,
            col_vals[i % ncol] + "o",
        )

        ax.plot(
            extdata[k].ext_waves["STIS"] * 1e3,
            extdata[k].ext_curve["STIS"] * norm_val + offset_val * i,
            col_vals[i % ncol] + "-",
        )

        if idlsave:
            (indxs,) = np.where(1.0 / extdata[k].ext_waves["MODEL"] < 9.0)
            ax.plot(
                extdata[k].ext_waves["MODEL"][indxs] * 1e3,
                extdata[k].ext_curve["MODEL"][indxs] * norm_val + offset_val * i,
                col_vals[i % ncol] + "--",
            )

    ax.set_xlim(90., 300.)
    ax.set_ylim(-1.0, 16.0 + offset_val * n_stars)

    ylimits = ax.get_ylim()
    xlimits = ax.get_xlim()
    ax.text(
        xlimits[0] + 0.25 * (xlimits[1] - xlimits[0]),
        ylimits[0] + 0.95 * (ylimits[1] - ylimits[0]),
        u"SMC, $E(B-V) > 0.15$",
        fontsize=1.5 * fontsize,
        horizontalalignment="left",
    )

    # ax.axvline(100., color='k')
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(ScalarFormatter())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--eps", help="save figure as an eps file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    fontsize = 16

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    fig, ax = pyplot.subplots(figsize=(8, 9))

#    datapath = "/home/kgordon/Hubble/SMCExt/Ed/"
    datapath = "/home/kgordon/Hubble/SMCExt/SENDTOKARL_FINAL/"
    plot_ext_stack(
        "smc_stars_reddened_good_highebv.dat",
        ax,
        idlsave=True,
        locpath=datapath,
        fontsize=fontsize,
    )
    ax.set_xlabel(r"$\lambda$ [$nm$]", fontsize=1.3 * fontsize)
    ax.set_ylabel(r"$E(\lambda - V)/E(B - V)$ + offset", fontsize=1.3 * fontsize)

    fig.tight_layout()

    if args.png:
        fig.savefig("smc_ext_goodonly.png")
    elif args.pdf:
        fig.savefig("smc_ext_goodonly.pdf")
    elif args.eps:
        fig.savefig("smc_ext_goodonly.eps")
    else:
        pyplot.show()
