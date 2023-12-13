import argparse

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable

from helpers import prettyname


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inclmc", help="include LMC points", action="store_true")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    fontsize = 12

    font = {"size": fontsize}

    plt.rc("font", **font)

    plt.rc("lines", linewidth=1)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("xtick.minor", width=2)
    plt.rc("ytick.major", width=2)
    plt.rc("ytick.minor", width=2)

    figsize = (8, 6)
    fig, ax = plt.subplots(figsize=figsize)

    tab1 = QTable.read("tables/qpah_from_image.dat", format="ascii.ipac")
    #tab1 = QTable.read("tables/samp_as_measured_ensemble_params.dat", format="ascii.ipac")
    tab2 = QTable.read("tables/gor24_smc_forecor_ensemble_params.dat", format="ascii.ipac")

    fname = ["highebv", "highebv_bumps"]  # , "highebv_flat"]  # , "lowebv"]
    fmts = ["bo", "rP"]  # , "cs"]  # , "vy"]
    legs = ["Weak/absent 2175 A bump", "Significant 2175 A bump"]  # , "Flat"]  # , r'$E(B-V)_\mathrm{SMC} < 0.1$']

    for cname, cfmt, cleg in zip(fname, fmts, legs):

        ptab = QTable.read(f"data/smc_stars_reddened_good_{cname}.dat", format="ascii", names=["name"], data_start=0)
        qpah = []
        qpah_unc = []
        B3 = []
        B3_unc = []
        for cname in ptab["name"]:
            mindx = np.where(prettyname(cname).replace(" ", "_") == tab1["name"])
            qpah.append(tab1["qpah"].data[mindx][0])
            qpah_unc.append(tab1["qpah_unc"].data[mindx][0])

            mindx = np.where(cname == tab2["name"])
            B3.append(tab2["B3"].data[mindx][0])
            B3_unc.append(tab2["B3_unc"].data[mindx][0])

        qpah = np.array(qpah)
        qpah_unc = np.array(qpah_unc)
        B3 = np.array(B3)
        B3_unc = np.array(B3_unc)
        ax.errorbar(qpah, B3, xerr=qpah_unc, yerr=B3_unc, fmt=cfmt, label=cleg)


    if args.inclmc:
        tab1 = QTable.read("tables/qpah_from_image_lmc.dat", format="ascii.ipac")
        tab2 = QTable.read("data/gor03_lmc_ensemble_params.dat", format="ascii.ipac")

        fname = ["avg", "lmc2"]  # , "highebv_flat"]  # , "lowebv"]
        fmts = ["gs", "mD"]  # , "cs"]  # , "vy"]
        legs = ["LMC Average Sample", "LMC2 Supershell Sample"]  # , "Flat"]  # , r'$E(B-V)_\mathrm{SMC} < 0.1$']

        for cname, cfmt, cleg in zip(fname, fmts, legs):

            ptab = QTable.read(f"data/lmc_positions_{cname}.dat", format="ascii.commented_header")
            qpah = []
            qpah_unc = []
            B3 = []
            B3_unc = []
            for cname in ptab["name"]:
                mindx = np.where(cname == tab1["name"])
                qpah.append(tab1["qpah"].data[mindx][0])
                qpah_unc.append(tab1["qpah_unc"].data[mindx][0])

                mindx = np.where(cname == tab2["Name"])
                C3 = tab2["C3"].data[mindx][0]
                C3_unc = tab2["C3_unc"].data[mindx][0]
                gamma = tab2["gamma"].value[mindx][0]
                gamma_unc = tab2["gamma_unc"].value[mindx][0]

                cB3 = C3 / gamma**2
                cB3_unc = cB3 * np.sqrt((C3_unc / C3)**2 + 2.*((gamma_unc / gamma)**2))
                B3.append(cB3)
                B3_unc.append(cB3_unc)

            qpah = np.array(qpah)
            qpah_unc = np.array(qpah_unc)
            B3 = np.array(B3)
            B3_unc = np.array(B3_unc)
            ax.errorbar(qpah, B3, xerr=qpah_unc, yerr=B3_unc, fmt=cfmt, label=cleg)

    ax.set_ylabel(r"B3 = 2175 $\mathrm{\AA}$ bump height")
    ax.set_xlabel(r"$q_\mathrm{PAH}$ = % dust mass in PAH grains")

    ax.legend()

    fig.tight_layout()

    save_str = "smcext_bump_vs_qpah"
    if args.png:
        fig.savefig(f"{save_str}.png")
    elif args.pdf:
        fig.savefig(f"{save_str}.pdf")
    else:
        plt.show()