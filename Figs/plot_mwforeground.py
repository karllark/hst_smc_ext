# plot the old and new curves for the same sightlines
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

from dust_extinction.parameter_averages import F19
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
    mpl.rc("lines", linewidth=1)
    mpl.rc("axes", linewidth=2)
    mpl.rc("xtick.major", width=2)
    mpl.rc("xtick.minor", width=2)
    mpl.rc("ytick.major", width=2)
    mpl.rc("ytick.minor", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), sharey=True, sharex=True)

    files = [
        "azv_018_ext_diff.fits",
        "azv_023_ext_diff.fits",
        "azv_214_ext_diff.fits",
        "azv_398_ext_diff.fits",
        "azv_456_ext_diff.fits",
    ]

    pcol = ["r", "b", "g", "c", "m"]
    allexts = []
    for k, cfile in enumerate(files):

        ext = ExtData(filename=f"fits/{cfile}")
        allexts.append(ext)

        ext.trans_elv_alav()
        print(cfile, ext.columns["AV"])
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

    ax[0].set_xlabel(r"$\lambda$ [$\mu$m$^{-1}$]")
    ax[1].set_xlabel(r"$\lambda$ [$\mu$m$^{-1}$]")
    ax[0].set_ylabel(r"$A(\lambda)/A(V)$")

    fig.tight_layout()

    # plot or save to a file
    outname = "mwforeground"
    if args.png:
        fig.savefig(outname + ".png")
    elif args.pdf:
        fig.savefig(outname + ".pdf")
    else:
        plt.show()
