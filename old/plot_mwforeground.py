# plot the old and new curves for the same sightlines
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import numpy as np

from dust_extinction.parameter_averages import F19, G16
from dust_extinction.averages import G03_SMCBar
from measure_extinction.extdata import ExtData, AverageExtData


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # plotting setup for easier to read plots
    fontsize = 18
    font = {"size": fontsize}
    mpl.rc("font", **font)
    mpl.rc("lines", linewidth=1.5)
    mpl.rc("axes", linewidth=2)
    mpl.rc("xtick.major", width=2)
    mpl.rc("xtick.minor", width=2)
    mpl.rc("ytick.major", width=2)
    mpl.rc("ytick.minor", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(14, 7), sharey=True, sharex=True)

    files = [
        "azv_018_ext_diff.fits",
        "azv_023_ext_diff.fits",
        "azv_214_ext_diff.fits",
        "azv_398_ext_diff.fits",
        "azv_456_ext_diff.fits",
    ]

    pcol = ["r", "b", "g", "c", "m"]
    allexts = []
    avs = []
    for k, cfile in enumerate(files):

        ext = ExtData(filename=f"fits/{cfile}")
        allexts.append(ext)

        ext.trans_elv_alav()
        avs.append(ext.columns["AV"])
        ext.plot(ax[0], wavenum=True, color=pcol[k], rebin_fac=None)

    ax[0].set_ylim(0.0, 7.0)

    # average the data
    aext = AverageExtData(allexts, min_number=3)
    aext.plot(ax[1], wavenum=True, color="b", rebin_fac=None)

    # plot MW expectations
    xmod = np.arange(0.3, 8.7, 0.1)
    fmod = F19(Rv=3.1)
    ax[1].plot(xmod, fmod(xmod), "k--")

    gmod = G03_SMCBar()
    ax[1].plot(xmod, gmod(xmod), "k:")

    gmod = G16(fA=0.6)
    ax[1].plot(xmod, gmod(xmod), "k-.")

    ax[0].set_xlabel(r"$\lambda$ [$\mu$m$^{-1}$]")
    ax[1].set_xlabel(r"$\lambda$ [$\mu$m$^{-1}$]")
    ax[0].set_ylabel(r"$A(\lambda)/A(V)$")

    snames = ["AzV18", "AzV23", "AzV214", "AzV398", "AzV456"]
    custom_lines = []
    for cname, ccol, cav in zip(snames, pcol, avs):
        custom_lines.append(Line2D([0], [0], color=ccol, lw=2, label=f"{cname} Av={cav:.2f}"))
    ax[0].legend(handles=custom_lines, fontsize=0.8 * fontsize)

    snames = ["Ave MW Foreground", "F19 Rv=3.1", "G03 SMCBar", "G16 fA=0.5"]
    pcol = ["b", "k", "k", "k"]
    pmark = ["-", "--", ":", "-."]
    custom_lines = []
    for cname, ccol, cmark in zip(snames, pcol, pmark):
        custom_lines.append(
            Line2D([0], [0], color=ccol, linestyle=cmark, lw=2, label=cname)
        )
    ax[1].legend(handles=custom_lines, fontsize=0.8 * fontsize)

    fig.tight_layout()

    # plot or save to a file
    outname = "mwforeground"
    if args.png:
        fig.savefig(outname + ".png")
    elif args.pdf:
        fig.savefig(outname + ".pdf")
    else:
        plt.show()
