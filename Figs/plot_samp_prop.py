import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib
from astropy.table import QTable

from measure_extinction.extdata import ExtData
from helpers import prettyname


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

    # gfilename = "data/smc_stars_reddened_good_highebv.dat"
    # lfilename = "data/smc_stars_reddened_good_lowebv.dat"
    # filename = "data/smc_stars_reddened_good.dat"
    filename = "data/smc_stars_all.dat"
    # filename = "data/smc_stars_reddened_suspect.dat"

    avs, avs_unc, ebvs, ebvs_unc, rvs, rvs_unc, nhs, nhs_unc, names = get_props(
        filename
    )

    # nhs from log to linear
    nhs_unc = 0.5 * (10 ** (nhs + nhs_unc) - 10 ** (nhs - nhs_unc))
    nhs = 10 ** nhs

    # get the MW HI foreground from radio measurements
    # mwfore = QTable.read("data/nhi_askap_karl.dat", format="ascii")
    mwfore = QTable.read("data/UVEXT_LMC_SMC_Coordinates_withHI_GASS.csv")

    # get the IR information
    irtab = QTable.read(
#        "data/forgordon_fit_param_smc.txt", format="ascii.commented_header"
        "data/forkarl_los_smc_fitprms_08mar23.txt", format="ascii.commented_header"
    )

    pnames = []
    nhs_forecor = np.copy(nhs)
    nhs_unc_forecor = np.copy(nhs_unc)
    ebvs_forecor = np.copy(ebvs)
    ebvs_unc_forecor = np.copy(ebvs_unc)
    avs_forecor = np.copy(avs)
    avs_unc_forecor = np.copy(avs_unc)
    nhs_mw = np.zeros(len(nhs))
    nhs_mw_unc = np.zeros(len(nhs))
    ebvs_mw = np.zeros(len(nhs))
    ebvs_mw_unc = np.zeros(len(nhs))
    qpahs = np.full((len(nhs)), np.NaN)
    for k, cname in enumerate(names):
        pname = prettyname(cname)
        pnames.append(pname)
        print(pname.replace(" ", "_"))
        (mindx,) = np.where(pname.replace(" ", "_") == mwfore["name"])
        if len(mindx) == 0:
            print(f"{cname} not found in MW foreground file")
        # nhs_mw[k] = 1e20 * mwfore["nhi_mw_askap"][mindx[0]]
        # nhs_mw_unc[k] = 1e20 * mwfore["nhi_mw_askap_err"][mindx[0]]
        nhs_mw[k] = 1e20 * mwfore["nhi"][mindx[0]]
        nhs_mw_unc[k] = 1e20 * mwfore["nhi_err"][mindx[0]]
        nhs_forecor[k] -= nhs_mw[k]
        nhs_unc_forecor[k] = np.sqrt(nhs_unc[k] ** 2 + nhs_mw_unc[k] ** 2)
        ebvs_mw[k] = nhs_mw[k] / 8.3e21
        ebvs_mw_unc[k] = nhs_mw_unc[k] / 8.3e21
        ebvs_forecor[k] -= ebvs_mw[k]
        avs_forecor[k] -= ebvs_mw[k] * 3.1
        avs_unc_forecor[k] = np.sqrt(avs_unc_forecor[k] ** 2 + (ebvs_mw_unc[k] * 3.1) ** 2)

        (mindx,) = np.where(pname.replace(" ", "") == irtab["Name"])
        if len(mindx) == 0:
            print(f"{cname} not found in IR properties file")
        else:
            qpahs[k] = irtab["q_PAH"][mindx[0]]  # in percentage terms

    # print out the names of low A(V) sightlines after MW foreground correction
    # names = np.array(names)
    # lnames = np.array(lnames)
    # print(names[avs_forecor < 0.2])
    # print(lnames[lavs_forecor < 0.2])

    print("E(B-V)_SMC", min(ebvs), max(ebvs), np.average(ebvs))
    print("E(B-V)_MW", min(ebvs_mw), max(ebvs_mw), np.average(ebvs_mw))

    # output a table of the sample properties
    outtab = QTable()
    outtab["name"] = pnames
    outtab["ebvs"] = ebvs
    outtab["ebvs_unc"] = ebvs_unc
    outtab["rvs"] = rvs
    outtab["rvs_unc"] = rvs_unc
    outtab["avs"] = avs
    outtab["avs_unc"] = avs_unc
    outtab["nhs"] = nhs
    outtab["nhs_unc"] = nhs_unc
    outtab["ebvs_forecor"] = ebvs_forecor
    outtab["ebvs_unc_forecor"] = ebvs_unc_forecor
    outtab["avs_forecor"] = avs_forecor
    outtab["avs_unc_forecor"] = avs_unc_forecor
    outtab["nhs_forecor"] = nhs_forecor
    outtab["nhs_unc_forecor"] = nhs_unc_forecor
    outtab["nhs_mw"] = nhs_mw
    outtab["nhs_mw_unc"] = nhs_mw_unc
    outtab["ebvs_mw"] = ebvs_mw
    outtab["ebvs_mw_unc"] = ebvs_mw_unc
    outtab["q_pahs"] = qpahs
    outtab.write(f"tables/samp_as_measured_ensemble_params.dat", format="ascii.ipac", overwrite=True)
    

    # outtab = outtab[np.argsort(outtab["name"])]
    outtab["nhs"] /= 1e20
    outtab["nhs_unc"] /= 1e20
    outtab["nhs_mw"] /= 1e20
    outtab["nhs_mw_unc"] /= 1e20
    outtab["nhs_forecor"] /= 1e20
    outtab["nhs_unc_forecor"] /= 1e20

    # loop over and write as a custom formatted latex table
    outfile = open("tables/derived_samp_properties.tex", "w")
    outfile2 = open("tables/ancillary_samp_properties.tex", "w")
    for row in outtab:
        outline = f"{row['name']}"
        #if row["avs_forecor"] / row["avs"] < 0.5:
        #    outline = f"{outline}*"
        #outline = f"{outline} & ${row['rvs']:.2f} \\pm {row['rvs_unc']:.2f}$"
        #outline = f"{outline} & ${row['avs']:.2f} \\pm {row['avs_unc']:.2f}$"
        outline = f"{outline} & ${row['nhs']:.2f} \\pm {row['nhs_unc']:.2f}$"
        outline = f"{outline} & ${row['ebvs']:.3f} \\pm {row['ebvs_unc']:.3f}$"
        outline = f"{outline} & ${row['nhs_mw']:.3f} \\pm {row['nhs_mw_unc']:.3f}$"
        outline = f"{outline} & ${row['ebvs_mw']:.3f} \\pm {row['ebvs_mw_unc']:.4f}$"
        #outline = f"{outline} & ${row['avs_forecor']:.2f} \\pm {row['avs_unc_forecor']:.2f}$"
        outline = f"{outline} & ${row['nhs_forecor']:.2f} \\pm {row['nhs_unc_forecor']:.2f}$"
        outline = f"{outline} & ${row['ebvs_forecor']:.3f} \\pm {row['ebvs_unc_forecor']:.3f}$"
        outfile.write(f"{outline} \\\\ \n")

        outline2 = f"{row['name']}"
        outline2 = f"{outline2} & ${row['nhs_mw']:.3f} \\pm {row['nhs_mw_unc']:.3f}$"
        outline2 = f"{outline2} & ${row['ebvs_mw']:.3f} \\pm {row['ebvs_mw_unc']:.4f}$"
        if row["q_pahs"] > 0.0:
            outline2 = f"{outline2} & ${row['q_pahs']:.2f}$"
        else:
            outline2 = f"{outline2} & \\nodata"
        outfile2.write(f"{outline2} \\\\ \n")
    outfile.close()
    outfile2.close()

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
        avs, rvs, xerr=avs_unc, yerr=rvs_unc, fmt="go", label="data",
    )
    ax[0].set_ylabel(r"$R(V)$")
    ax[0].set_xlabel(r"$A(V)$")
    ax[0].tick_params("both", length=10, width=2, which="major")
    ax[0].tick_params("both", length=5, width=1, which="minor")

    # N(HI) versus E(B-V)
    ax[1].errorbar(
        ebvs,
        nhs,
        xerr=ebvs_unc,
        yerr=nhs_unc,
        fmt="go",
        markerfacecolor="none",
        alpha=0.5,
        # label="E(B-V) > 0.15",
    )
    ax[1].errorbar(
        ebvs_forecor,
        nhs_forecor,
        xerr=ebvs_unc_forecor,
        yerr=nhs_unc,
        fmt="go",
        alpha=0.8,
    #     # label="E(B-V) > 0.15",
    )
    ax[1].set_xlabel(r"$E(B-V)$")
    ax[1].set_ylabel(r"$N(HI)$")
    ax[1].tick_params("both", length=10, width=2, which="major")
    ax[1].tick_params("both", length=5, width=1, which="minor")

    # plot MW and SMC expected lines
    ax[1].plot(
        [0.0, .5], [0.0, 13.18e21 * 3.1 * .5], "k--", label="SMC (Gordon et al. 2003)"
    )
    ax[1].plot([0.0, .5], [0.0, 1.55e21 * 2.74 * .5], "k:", label="MW (Bohlin et al. 1978+)")

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
