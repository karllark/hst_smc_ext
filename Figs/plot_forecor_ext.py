import argparse
import copy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import astropy.units as u

from measure_extinction.extdata import ExtData
from dust_extinction.parameter_averages import G23


def foreground_correct_extinction(ext, forehi, forehi_unc, foremod):
    """
    Correct an extinction curve for foreground extinction
    """
    # foreground corrected extinction
    ext_fc = copy.deepcopy(ext)
    foreebv = forehi / 5e21
    foreebv_unc = foreebv * (forehi_unc / forehi)
    foreax = foreebv * 3.1
    foreax_unc = foreax * (forehi_unc / forehi)

    for src in ext.waves.keys():
        # get foreground ext in A(l)/A(V)
        foreext = foremod(ext.waves[src])
        if ext.type == "elx":
            foreext_elv = (foreext - 1.0) * foreax
            ext_fc.exts[src] -= foreext_elv
            ext_fc.columns["AV"] = (
                ext.columns["AV"][0] - foreax,
                np.sqrt(ext_fc.columns["AV"][1] ** 2 + foreax_unc ** 2),
            )
            ext_fc.columns["EBV"] = (
                ext.columns["EBV"][0] - foreebv,
                np.sqrt(ext_fc.columns["EBV"][1] ** 2 + foreebv_unc ** 2),
            )
            rv = ext_fc.columns["AV"][0] / ext_fc.columns["EBV"][0]
            rv_unc = ((ext_fc.columns["AV"][1] / ext_fc.columns["AV"][0]) ** 2
                      + (ext_fc.columns["EBV"][1] / ext_fc.columns["EBV"][0]) ** 2)
            rv_unc = rv * np.sqrt(rv_unc)
            ext_fc.columns["RV"] = (rv, rv_unc)
            if ext.columns["LOGHI"][0] > 0.0:
                hi = 10 ** ext.columns["LOGHI"][0] - forehi
                hi_up = 10 ** (ext.columns["LOGHI"][0] + ext.columns["LOGHI"][1]) + forehi_unc
                hi_down = 10 ** (ext.columns["LOGHI"][0] - ext.columns["LOGHI"][1]) - forehi_unc
                if hi > 0:
                    loghi = np.log10(hi)
                    loghi_unc = (np.log10(hi_up) - np.log10(hi_down)) / 2.
                else:
                    loghi = 0.0
                    loghi_unc = 0.0
            else:
                loghi = 0.0
                loghi_unc = 0.0

            ext_fc.columns["LOGHI"] = (loghi, loghi_unc)
        else:
            print(f"{ext.type} not supported.")
            exit()

    return ext_fc


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("extname", help="extinction curve name")
    parser.add_argument("forehi", help="foreground HI column [10^21]", type=float)
    parser.add_argument("forehi_unc", help="foreground HI column unc [10^21]", type=float)
    parser.add_argument("--adjusted", help="foreground HI column has been adjusted", action="store_true")
    parser.add_argument("--prev", help="previous extinction")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    args.forehi *= 1e21
    args.forehi_unc *= 1e21

    if "mr12" in args.extname:
        rebinfac = 10
    else:
        rebinfac = None

    fontsize = 14

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    figsize = (8.0, 8.0)
    fig, ax = plt.subplots(nrows=2, figsize=figsize)

    ext = ExtData(filename=f"fits/{args.extname}_ext.fits")
    #tot_av = ext.columns["EBV"][0] * ext.columns["RV"][0]
    #ext.columns["AV"] = (tot_av, 0.0)

    # do foreground correction
    extmod = G23(Rv=3.1)
    ext_fc = foreground_correct_extinction(ext, args.forehi, args.forehi_unc, extmod)
    savename = f"fits/{args.extname}_ext_forecor.fits"
    if args.adjusted:
        savename = savename.replace(".fits", "_adjusted.fits")
    ext_fc.save(savename)

    # make plots
    ext.plot(
        ax[0],
        color="b",
        legend_key="IUE",
        legend_label="SMC+MW foreground",
        rebin_fac=rebinfac,
        wavenum=True,
    )

    ext_fc.plot(
        ax[0],
        color="g",
        legend_key="IUE",
        legend_label="SMC only",
        rebin_fac=rebinfac,
        wavenum=True,
    )

    ext.trans_elv_elvebv()
    ext.plot(
        ax[1],
        color="b",
        legend_key="IUE",
        legend_label="SMC+MW foreground",
        rebin_fac=rebinfac,
        wavenum=True,
    )

    ext_fc.trans_elv_elvebv()
    ext_fc.plot(
        ax[1],
        color="g",
        legend_key="IUE",
        legend_label="SMC only",
        rebin_fac=rebinfac,
        wavenum=True,
    )

    if args.prev:
        pext = ExtData(filename=args.prev)
        if "IUE" in pext.waves.keys():
            legend_key = "IUE"
        else:
            legend_key = "STIS"

        pext.plot(ax[0], color="m", legend_key=legend_key, legend_label="G03", alpha=0.5, wavenum=True)
        pext.trans_elv_elvebv()
        pext.plot(ax[1], color="m", legend_key=legend_key, legend_label="G03", alpha=0.5, wavenum=True)

    # add in average MW with R(V) = 3.1
    mwaves = np.logspace(np.log10(0.1), np.log10(3.0), num=500) * u.micron
    ax[1].plot(1.0 / mwaves, (extmod(mwaves) - 1.0) * 3.1, "k-", label="G23 R(V)=3.1")

    ax[0].set_xlabel("")
    ax[0].tick_params("both", length=10, width=2, which="major")
    ax[0].tick_params("both", length=5, width=1, which="minor")

    ax[0].set_xlim(0.0, 9.0)
    ax[1].set_xlim(0.0, 9.0)
    ax[0].set_ylim(-1.0, 4.0)
    ax[0].set_ylabel(r"E($\lambda$ - V)")

    ax[1].set_ylim(-4.0, 18.0)
    #ax[1].set_xscale("log")

    tot_ebv = ext.columns["EBV"][0]
    ax[0].legend(title=f"{args.extname}\ntotal E(B-V)={tot_ebv}\nforeground E(B-V)={args.forehi / 5e21:.3f}",
                 fontsize=0.8*fontsize)
    ax[1].legend(title=rf"{args.extname}", fontsize=0.8*fontsize)

    # for 2nd x-axis with lambda values
    axis_xs = np.array([0.12, 0.15, 0.2, 0.3, 0.5, 1.0, 2.0])
    new_ticks = 1 / axis_xs
    new_ticks_labels = ["%.2f" % z for z in axis_xs]
    tax = ax[0].twiny()
    tax.set_xlim(ax[0].get_xlim())
    tax.set_xticks(new_ticks)
    tax.set_xticklabels(new_ticks_labels, fontsize=0.8*fontsize)
    tax.set_xlabel(r"$\lambda$ [$\mu$m]")

    fig.tight_layout()

    save_str = f"fits/forecor_ext_{args.extname}"
    if args.adjusted:
        save_str = f"{save_str}_adjusted"
    if args.png:
        fig.savefig(f"{save_str}.png")
    elif args.pdf:
        fig.savefig(f"{save_str}.pdf")
    else:
        plt.show()
