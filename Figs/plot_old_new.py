# plot the old and new curves for the same sightlines
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from measure_extinction.extdata import ExtData


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

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

    fig, fax = plt.subplots(nrows=2, ncols=3, figsize=(14, 10), sharey=True, sharex=True)

    ax = [fax[0, 0], fax[0, 1], fax[1, 0], fax[1, 1], fax[0, 2]]

    ptypes = ["IUE", "BAND"]
    psym = ["-", "o"]

    snames = ["AzV18", "AzV23", "AzV214", "AzV398", "AzV456"]
    ofiles = [
        "azv18_azv462_ext.fits",
        "azv23_azv404_ext.fits",
        "azv214_azv380_ext.fits",
        "azv398_azv289_ext.fits",
        "azv456_azv70_ext.fits",
    ]
    nfiles = [
        "azv_018_ext.fits",
        "azv_023_ext.fits",
        "azv_214_ext.fits",
        "azv_398_ext.fits",
        "azv_456_ext.fits",
    ]
    for k, ofile in enumerate(ofiles):

        oext = ExtData(filename=f"prev/{ofile}")
        next = ExtData(filename=f"fits/{nfiles[k]}")

        for curtype, csym in zip(ptypes, psym):
            oext.exts[curtype][oext.npts[curtype] == 0] = np.nan
            if curtype == "IUE":
                olabel = "G03"
                nlabel = "this work"
                dlabel = "this work - G03"
                dlabel2 = f"{snames[k]} diff"
            else:
                olabel = None
                nlabel = None
                dlabel = None
                dlabel2 = None
            ax[k].plot(
                1.0 / oext.waves[curtype], oext.exts[curtype], f"g{csym}", label=olabel
            )
            next.exts[curtype][next.npts[curtype] == 0] = np.nan
            ax[k].plot(
                1.0 / next.waves[curtype],
                next.exts[curtype],
                f"b{csym}",
                label=nlabel,
            )

            # plot the difference between new-old -> foreground
            iold = np.interp(
                next.waves[curtype], oext.waves[curtype], oext.exts[curtype]
            )
            ax[k].plot(
                1.0 / next.waves[curtype],
                next.exts[curtype] - iold,
                f"r{csym}",
                label=dlabel,
            )
            fax[1, 2].plot(
                1.0 / next.waves[curtype],
                next.exts[curtype] - iold,
                label=dlabel2,
            )

        ax[k].set_title(snames[k])
        # ax[k].legend()

    ax[0].legend()
    fax[1, 2].legend()

    ax[0].set_ylim(-2.0, 5.0)

    fax[1, 0].set_xlabel(r"$\lambda$ [$\mu$m$^{-1}$]")
    fax[1, 1].set_xlabel(r"$\lambda$ [$\mu$m$^{-1}$]")
    fax[1, 2].set_xlabel(r"$\lambda$ [$\mu$m$^{-1}$]")

    fax[0, 0].set_ylabel(r"$E(\lambda - V)$")
    fax[1, 0].set_ylabel(r"$E(\lambda - V)$")

    fig.tight_layout()

    # plot or save to a file
    outname = "old_vs_new"
    if args.png:
        fig.savefig(outname + ".png")
    elif args.pdf:
        fig.savefig(outname + ".pdf")
    else:
        plt.show()
