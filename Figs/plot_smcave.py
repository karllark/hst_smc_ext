import argparse
import copy

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
import astropy.units as u

from dust_extinction.averages import G03_SMCBar
from dust_extinction.parameter_averages import G23

from measure_extinction.merge_obsspec import merge_stis_obsspec
from measure_extinction.extdata import ExtData, AverageExtData


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bumps", help="save figure as a png file", action="store_true")    
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

    filelist = "data/smc_stars_reddened_good_highebv.dat"
    if args.bumps:
        filelist = filelist.replace(".dat", "_bumps.dat")
        klabel = "SMC Bumps"
        ofilename = "fits/smcave_bumps.fits"
    else:
        klabel = "SMC Average"
        ofilename = "fits/smcave.fits"
    f = open(filelist, "r")
    file_lines = list(f)
    starnames = []
    extdatas = []
    spslopes = []
    rebinfac = None
    rebinfac = 5
    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            starnames.append(name)
            text = ExtData()
            filebase = "./fits/" + starnames[-1] + "_ext_forecor"
            filebase = f"{filebase}_FM90"
            text.read(f"{filebase}.fits")

            # put all the UV spectra onto the same wavelength grid
            #   mock up the data as an observed STIS spectrum and use merge_stis_obsspec
            if "IUE" in text.waves.keys():
                uvtype = "IUE"
            else:
                uvtype = "STIS"
            ctable = QTable()
            ctable["WAVELENGTH"] = np.copy(text.waves[uvtype].to(u.Angstrom))
            ctable["FLUX"] = np.copy(text.exts[uvtype])
            ctable["STAT-ERROR"] = np.copy(text.uncs[uvtype])
            ctable["NPTS"] = np.full(len(text.waves[uvtype]), 1.)

            ntable = merge_stis_obsspec([ctable])
            next = copy.deepcopy(text)
            if "IUE" in next.waves.keys():
                del next.waves["IUE"]  # remove so it is not averaged/plotted
                del next.exts["IUE"]
                del next.uncs["IUE"]
                del next.npts["IUE"]
            next.waves["STIS"] = ntable["WAVELENGTH"]
            next.exts["STIS"] = ntable["FLUX"]
            next.uncs["STIS"] = ntable["SIGMA"]
            next.npts["STIS"] = ntable["NPTS"]

            # get R(V) and transform to A(lambda)/A(V)
            # rv = text.columns["RV"][0]
            # next.exts["BAND"] = next.exts["BAND"] / rv + 1
            # next.exts["STIS"] = next.exts["STIS"] / rv + 1

            extdatas.append(next)

            next.plot(
                ax,
                color="k",
                alpha=0.1,
                rebin_fac=rebinfac,
                wavenum=True,
            )
    aveext = AverageExtData(extdatas, min_number=2)

    # calculates R(V) as E(lambda-V)/E(B-V) is input
    aveext.calc_AV_JHK()
    rv = aveext.columns["AV"][0]
    print(f"Average has R(V) = {rv:.2f}")
    #aveext.exts["BAND"] = aveext.exts["BAND"] / rv + 1
    #aveext.exts["STIS"] = aveext.exts["STIS"] / rv + 1
    aveext.plot(ax, color="b", rebin_fac=rebinfac, wavenum=True, legend_key="STIS", legend_label=klabel)

    ax.plot(G03_SMCBar.obsdata_x, (G03_SMCBar.obsdata_axav - 1)*G03_SMCBar.Rv, 'ko',
        label="G03 SMCBar")
    aveext.save(ofilename)

    mod_x = np.arange(0.5, 9.0, 0.1)
    mwmod = G23()
    mod_y = (mwmod(mod_x) - 1) * 3.1
    ax.plot(mod_x, mod_y, "g--", alpha=0.75, label="Milky Way: G23 R(V)=3.1", linewidth=3.0)

    ax.set_ylim(-5.0, 20.0)
    ax.legend()

    fig.tight_layout()

    save_str = "smcext_ave"
    if args.bumps:
        save_str = f"{save_str}_bumps"
    if args.png:
        fig.savefig(f"{save_str}.png")
    elif args.pdf:
        fig.savefig(f"{save_str}.pdf")
    else:
        plt.show()