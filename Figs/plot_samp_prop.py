import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib
from astropy.table import QTable

from measure_extinction.extdata import ExtData


def get_props(filename):

    f = open(filename, "r")
    file_lines = list(f)
    extnames = []
    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            extnames.append(name)

    n_ext = len(extnames)
    avs = np.full((n_ext), 0.0)
    avs_unc = np.full((n_ext), 0.0)
    ebvs = np.full((n_ext), 0.0)
    ebvs_unc = np.full((n_ext), 0.0)
    rvs = np.full((n_ext), 0.0)
    rvs_unc = np.full((n_ext), 0.0)
    nhs = np.full((n_ext), 0.0)
    nhs_unc = np.full((n_ext), 0.0)

    for k, cname in enumerate(extnames):

        bfile = f"fits/{cname}_ext.fits"
        cext = ExtData(filename=bfile)

        ebv = cext.columns["EBV"][0]
        ebv_unc = cext.columns["EBV"][1]
        rv = cext.columns["RV"][0]
        rv_unc = cext.columns["RV"][1]
        nh = cext.columns["LOGHI"][0]
        nh_unc = cext.columns["LOGHI"][1]
        av = rv * ebv
        av_unc = av * np.sqrt((rv_unc / rv) ** 2 + (ebv_unc / ebv) ** 2)

        ebvs[k] = ebv
        ebvs_unc[k] = ebv_unc
        rvs[k] = rv
        rvs_unc[k] = rv_unc
        avs[k] = av
        avs_unc[k] = av_unc
        nhs[k] = nh
        nhs_unc[k] = nh_unc

    return (avs, avs_unc, ebvs, ebvs_unc, rvs, rvs_unc, nhs, nhs_unc, extnames)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    gfilename = "data/smc_stars_reddened_good_highebv.dat"
    lfilename = "data/smc_stars_reddened_good_lowebv.dat"

    avs, avs_unc, ebvs, ebvs_unc, rvs, rvs_unc, nhs, nhs_unc, names = get_props(
        gfilename
    )
    (
        lavs,
        lavs_unc,
        lebvs,
        lebvs_unc,
        lrvs,
        lrvs_unc,
        lnhs,
        lnhs_unc,
        lnames,
    ) = get_props(lfilename)

    # nhs from log to linear
    nhs_unc = 0.5 * (10 ** (nhs + nhs_unc) - 10 ** (nhs - nhs_unc))
    nhs = 10 ** nhs
    lnhs_unc = 0.5 * (10 ** (lnhs + lnhs_unc) - 10 ** (lnhs - lnhs_unc))
    lnhs = 10 ** lnhs

    # get the MW HI foreground from radio measurements
    mwfore = QTable.read("data/nhi_askap_karl.dat", format="ascii")

    nhs_forecor = np.copy(nhs)
    avs_forecor = np.copy(avs)
    for k, cname in enumerate(names):
        (mindx,) = np.where(cname == mwfore["name"])
        if len(mindx) == 0:
            print(f"{cname} not found in MW foreground file")
        mwforeground = 1e20 * mwfore["nhi_mw_askap"][mindx[0]]
        nhs_forecor[k] -= mwforeground
        avs_forecor[k] -= mwforeground / 1.55e21

    lnhs_forecor = np.copy(lnhs)
    lavs_forecor = np.copy(lavs)
    for k, cname in enumerate(lnames):
        (mindx,) = np.where(cname == mwfore["name"])
        if len(mindx) == 0:
            print(f"{cname} not found in MW foreground file")
        mwforeground = 1e20 * mwfore["nhi_mw_askap"][mindx[0]]
        lnhs_forecor[k] -= mwforeground
        lavs_forecor[k] -= mwforeground / 1.55e21

    # print out the names of low A(V) sightlines after MW foreground correction
    names = np.array(names)
    lnames = np.array(lnames)
    print(names[avs_forecor < 0.2])
    print(lnames[lavs_forecor < 0.2])

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

    figsize = (10.0, 5.0)
    fig, ax = pyplot.subplots(ncols=2, figsize=figsize)

    # R(V) versus A(V)
    ax[0].errorbar(
        avs, rvs, xerr=avs_unc, yerr=rvs_unc, fmt="go", label="E(B-V) > 0.15",
    )
    ax[0].errorbar(
        lavs, lrvs, xerr=lavs_unc, yerr=lrvs_unc, fmt="bo", label="E(B-V) < 0.15",
    )
    ax[0].set_ylabel(r"$R(V)$")
    ax[0].set_xlabel(r"$A(V)$")
    ax[0].tick_params("both", length=10, width=2, which="major")
    ax[0].tick_params("both", length=5, width=1, which="minor")

    # N(HI) versus A(V)
    ax[1].errorbar(
        avs,
        nhs,
        xerr=avs_unc,
        yerr=nhs_unc,
        fmt="go",
        # label="E(B-V) > 0.15",
    )
    ax[1].errorbar(
        avs_forecor,
        nhs_forecor,
        xerr=avs_unc,
        yerr=nhs_unc,
        fmt="go",
        alpha=0.3,
        # label="E(B-V) > 0.15",
    )
    ax[1].errorbar(
        lavs,
        lnhs,
        xerr=lavs_unc,
        yerr=lnhs_unc,
        fmt="bo",
        # label="E(B-V) < 0.15",
    )
    ax[1].errorbar(
        lavs_forecor,
        lnhs_forecor,
        xerr=lavs_unc,
        yerr=lnhs_unc,
        fmt="bo",
        alpha=0.3,
        # label="E(B-V) < 0.15",
    )
    ax[1].set_xlabel(r"$A(V)$")
    ax[1].set_ylabel(r"$N(HI)$")
    ax[1].tick_params("both", length=10, width=2, which="major")
    ax[1].tick_params("both", length=5, width=1, which="minor")

    # plot MW and SMC expected lines
    ax[1].plot(
        [0.0, 1.5], [0.0, 13.18e21 * 1.5], "k--", label="SMC (Gordon et al. 2003)"
    )
    ax[1].plot([0.0, 1.5], [0.0, 1.55e21 * 1.5], "k:", label="MW (Bohlin et al. 1978+)")

    # legends
    ax[0].legend(fontsize=0.8 * fontsize)
    ax[1].legend(fontsize=0.8 * fontsize)

    fig.tight_layout()

    save_str = "smcext_sampprop"
    if args.png:
        fig.savefig(f"{save_str}.png")
    elif args.pdf:
        fig.savefig(f"{save_str}.pdf")
    else:
        pyplot.show()
