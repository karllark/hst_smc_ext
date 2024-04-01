import argparse
import math

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
from astropy.table import QTable
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip

from helpers import prettyname
from fit_full2dcor import lnlike_correlated


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inclmc", help="include LMC points", action="store_true")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    fontsize = 16

    font = {"size": fontsize}

    plt.rc("font", **font)

    plt.rc("lines", linewidth=1)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("xtick.minor", width=2)
    plt.rc("ytick.major", width=2)
    plt.rc("ytick.minor", width=2)

    figsize = (14, 6)
    fig, fax = plt.subplots(ncols=2, figsize=figsize)
    ax = fax[0]
    ax2 = fax[1]

    tab1 = QTable.read("tables/qpah_from_image.dat", format="ascii.ipac")
    #tab1 = QTable.read("tables/samp_as_measured_ensemble_params.dat", format="ascii.ipac")
    tab2 = QTable.read("tables/gor24_smc_ensemble_params.dat", format="ascii.ipac")

    fname = ["highebv", "highebv_bumps"]  # , "highebv_flat"]  # , "lowebv"]
    fmts = ["bo", "rP"]  # , "cs"]  # , "vy"]
    legs = ["Weak/absent 2175 A bump", "Significant 2175 A bump"]  # , "Flat"]  # , r'$E(B-V)_\mathrm{SMC} < 0.1$']

    xvals = None
    for cname, cfmt, cleg in zip(fname, fmts, legs):

        ptab = QTable.read(f"data/smc_stars_reddened_good_{cname}.dat", format="ascii", names=["name"], data_start=0)
        qpah = []
        qpah_unc = []
        B3 = []
        B3_unc = []
        C4 = []
        C4_unc = []
        for cname in ptab["name"]:
            mindx = np.where(prettyname(cname).replace(" ", "_") == tab1["name"])
            qpah.append(tab1["qpah"].data[mindx][0])
            qpah_unc.append(tab1["qpah_unc"].data[mindx][0])

            mindx = np.where(cname == tab2["name"])
            area = 0.5 * math.pi * tab2["B3"].data[mindx][0] * tab2["gamma"].data[mindx][0]
            area_unc = area * np.sqrt((tab2["B3_unc"].data[mindx][0] / tab2["B3"].data[mindx][0])**2
                                      + (tab2["gamma_unc"].data[mindx][0] / tab2["gamma"].data[mindx][0])**2)
            B3.append(area)
            B3_unc.append(np.absolute(area_unc))

            C4.append(tab2["C4"].data[mindx][0])
            C4_unc.append(tab2["C4_unc"].data[mindx][0])

        qpah = np.array(qpah)
        qpah_unc = np.array(qpah_unc)
        B3 = np.array(B3)
        B3_unc = np.array(B3_unc)
        ax.errorbar(qpah, B3, xerr=qpah_unc, yerr=B3_unc, fmt=cfmt, label=cleg)
        ax2.errorbar(qpah, C4, xerr=qpah_unc, yerr=C4_unc, fmt=cfmt, label=cleg)

        if xvals is not None:
            xvals = np.concatenate((xvals, qpah))
            xvals_unc = np.concatenate((xvals_unc, qpah_unc))
            yvals = np.concatenate((yvals, B3))
            yvals_unc = np.concatenate((yvals_unc, B3_unc))
            yvals2 = np.concatenate((yvals2, C4))
            yvals2_unc = np.concatenate((yvals2_unc, C4_unc))
        else:
            xvals = qpah
            xvals_unc = qpah_unc
            yvals = B3
            yvals_unc = B3_unc
            yvals2 = C4
            yvals2_unc = C4_unc


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
            C4 = []
            C4_unc = []
            for cname in ptab["name"]:
                mindx = np.where(cname == tab1["name"])
                qpah.append(tab1["qpah"].data[mindx][0])
                qpah_unc.append(tab1["qpah_unc"].data[mindx][0])

                mindx = np.where(cname == tab2["Name"])
                C3 = tab2["C3"].data[mindx][0]
                C3_unc = tab2["C3_unc"].data[mindx][0]
                gamma = tab2["gamma"].value[mindx][0]
                gamma_unc = tab2["gamma_unc"].value[mindx][0]

                area = 0.5 * math.pi * C3 / gamma
                area_unc = area * np.sqrt((C3_unc / C3)**2 + 2.*((gamma_unc / gamma)**2))
                B3.append(area)
                B3_unc.append(area_unc)

                C4.append(tab2["C4"].data[mindx][0])
                C4_unc.append(tab2["C4_unc"].data[mindx][0])

            qpah = np.array(qpah)
            qpah_unc = np.array(qpah_unc)
            B3 = np.array(B3)
            B3_unc = np.array(B3_unc)
            ax.errorbar(qpah, B3, xerr=qpah_unc, yerr=B3_unc, fmt=cfmt, label=cleg, alpha=0.4)
            ax2.errorbar(qpah, C4, xerr=qpah_unc, yerr=C4_unc, fmt=cfmt, label=cleg, alpha=0.4)

            xvals = np.concatenate((xvals, qpah))
            xvals_unc = np.concatenate((xvals_unc, qpah_unc))
            yvals = np.concatenate((yvals, B3))
            yvals_unc = np.concatenate((yvals_unc, B3_unc))
            yvals2 = np.concatenate((yvals2, C4))
            yvals2_unc = np.concatenate((yvals2_unc, C4_unc))

    ax.set_ylabel(r"$\pi B_3 \gamma$ / 2 = 2175 $\mathrm{\AA}$ bump area")
    ax.set_xlabel(r"$q_\mathrm{PAH}$ = % dust mass in PAH grains")

    #ax2.yaxis.tick_right()
    #ax2.yaxis.set_label_position("right")
    ax2.set_ylabel(r"C4 = FUV rise amplitude")
    ax2.set_xlabel(r"$q_\mathrm{PAH}$ = % dust mass in PAH grains")

    fitline = True
    if fitline:

        # now fit a line to the data
        npts = len(xvals)
        covs = np.zeros((npts, 2, 2))
        covs2 = np.zeros((npts, 2, 2))
        for k in range(npts):
            covs[k, 0, 0] = xvals_unc[k]** 2
            covs[k, 0, 1] = 0.0
            covs[k, 1, 0] = 0.0
            covs[k, 1, 1] = yvals_unc[k]**2

            covs2[k, 0, 0] = xvals_unc[k]**2
            covs2[k, 0, 1] = 0.0
            covs2[k, 1, 0] = 0.0
            covs2[k, 1, 1] = yvals2_unc[k]**2
            if not np.all(np.linalg.eigvals(covs[k, :, :]) > 0):
                print(k, np.all(np.linalg.eigvals(covs[k, :, :]) > 0))
                print(covs[k, :, :])
        intinfo = [-5.0, 15.0, 0.1]

        def nll(*args):
            return -lnlike_correlated(*args)

        line_orig = models.Linear1D(slope=1.0, intercept=0.0)
        fit = fitting.LinearLSQFitter()
        or_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip, niter=3, sigma=3.0)
        x = np.arange(0.0, 7., 0.1)

        # fit B3 area
        fitted_line, mask = or_fit(line_orig, xvals, yvals, weights=1/yvals_unc)

        masked_data = np.ma.masked_array(yvals, mask=~mask)
        ax.plot(xvals, masked_data, "ko", fillstyle="none", ms=10, label="Not used in fit")

        result = op.minimize(nll, fitted_line.parameters, args=(yvals[~mask], fitted_line, covs, intinfo, xvals[~mask]))
        nparams = result["x"]
        print("B3 area fit: ", nparams)

        fitted_line2 = models.Linear1D(slope=nparams[0], intercept=nparams[1])
        ax.plot(x, fitted_line2(x), "k--", label="Linear Fit")

        # fit C4
        fitted_line, mask = or_fit(line_orig, xvals, yvals2, weights=1/yvals2_unc)

        masked_data = np.ma.masked_array(yvals2, mask=~mask)
        ax2.plot(xvals, masked_data, "ko", fillstyle="none", ms=10, label="Not used in fit")

        result = op.minimize(nll, fitted_line.parameters, args=(yvals2[~mask], fitted_line, covs, intinfo, xvals[~mask]))
        nparams = result["x"]
        print("C4 fit: ", nparams)

        fitted_line2 = models.Linear1D(slope=nparams[0], intercept=nparams[1])
        ax2.plot(x, fitted_line2(x), "k--", label="Linear Fit")

        #get handles and labels
        handles, labels = ax.get_legend_handles_labels()

        #specify order of items in legend
        order = [2, 3, 4, 5, 1, 0]

        #add legend to plot
        ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize=0.8*fontsize)
    else:
        ax.legend(fontsize=0.8*fontsize)

    fig.tight_layout()

    save_str = "smcext_bump_vs_qpah"
    if args.png:
        fig.savefig(f"{save_str}.png")
    elif args.pdf:
        fig.savefig(f"{save_str}.pdf")
    else:
        plt.show()