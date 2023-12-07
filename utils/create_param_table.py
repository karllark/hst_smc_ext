import glob

import numpy as np
from astropy.table import QTable

# import astropy.units as u

from measure_extinction.extdata import ExtData

from helpers import prettyname


if __name__ == "__main__":

    for ctype in ["", "_forecor"]:

        otab = QTable(
            # fmt: off
            names=("name", "AV", "AV_unc", "EBV", "EBV_unc", "RV", "RV_unc", "LOGHI", "LOGHI_unc",
                   "C1", "C1_unc", "C2", "C2_unc", "B3", "B3_unc", "C4", "C4_unc",
                   "x0", "x0_unc", "gamma", "gamma_unc"),
            dtype=("S", "f", "f", "f", "f", "f", "f", "f", "f", "f", "f", "f",
                   "f", "f", "f", "f", "f", "f", "f", "f", "f"),
            # fmt:on
        )

        otab_lat = QTable(
            # fmt: off
            names=("Name", 
                   "$E(B-V)$", "$R(V)$", "$N_\mathrm{SMC}(HI)$",
                   "$C_1$", "$C_2$", "$B_3$", "$C_4$"),
            dtype=("S", "S", "S", "S", "S", "S", "S", "S")
            # fmt:on
        )

        colnames = ["AV", "EBV", "RV", "LOGHI"]
        fm90names = ["C1", "C2", "B3", "C4", "XO", "GAMMA"]

        files = ["highebv_bumps", "highebv", "lowebv"]
        tags = [r"$E(B-V)_\mathrm{SMC} \geq 0.1, Significant Bump",
                r"$E(B-V)_\mathrm{SMC} \geq 0.1, Weak/Absent Bump",
                r"$E(B-V)_\mathrm{SMC} < 0.1, Weak/Absent Bump"]
        for cfile, ctag in zip(files, tags):

            filename = f"data/smc_stars_reddened_good_{cfile}.dat"

            f = open(filename, "r")
            file_lines = list(f)
            starnames = []
            for line in file_lines:
                if (line.find("#") != 0) & (len(line) > 0):
                    name = line.rstrip()
                    starnames.append(name)
            # starnames = np.sort(starnames)

            for cname in starnames:
                print(cname)
                cfile = f"fits/{cname}_ext{ctype}_FM90.fits"

                edata = ExtData(filename=cfile)
                if "AV" not in edata.columns:
                    av = edata.columns["EBV"][0] * edata.columns["RV"][0]
                    # fmt: off
                    av_unc = ((edata.columns["EBV"][0] / edata.columns["EBV"][0]) ** 2 
                            + (edata.columns["RV"][0] / edata.columns["RV"][0]) ** 2)
                    # fmt: on
                    av_unc = av * np.sqrt(av_unc)
                    edata.columns["AV"] = (av, av_unc)

                # get foreground subtraction plus / minus fits
                if ctype == "_forecor":
                    pedata = ExtData(filename=cfile.replace("forecor", "forecor_plus"))
                    medata = ExtData(filename=cfile.replace("forecor", "forecor_minus"))
                    dav = 0.5 * abs(medata.columns["AV"][0] - pedata.columns["AV"][0])
                    # print(cname, edata.columns["AV"], dav)

                rdata = []
                rdata_lat = []
                rdata.append(cname)
                pcname = prettyname(cname)
                rdata_lat.append(pcname)
                for ccol in colnames:
                    val = edata.columns[ccol][0]
                    unc = edata.columns[ccol][1]
                    if ctype == "_forecor":
                        # fmt: off
                        unc = np.sqrt((unc ** 2) 
                                    + (np.absolute(medata.columns[ccol][0] - pedata.columns[ccol][0])) ** 2)
                        # fmt: on
                    rdata.append(val)
                    rdata.append(unc)
                    if ccol != "AV":
                        rdata_lat.append(fr"${val:.2f} \pm {unc:.2f}$")

                for ccol in fm90names:
                    val = edata.fm90_p50_fit[ccol][0]
                    unc = edata.fm90_p50_fit[ccol][1]
                    if ctype == "_forecor":
                        # fmt: off
                        unc = np.sqrt((unc ** 2) 
                                    + (np.absolute(medata.fm90_p50_fit[ccol][0] - pedata.fm90_p50_fit[ccol][0])) ** 2)
                        # fmt: on
                    rdata.append(val)
                    rdata.append(unc)
                    if ccol not in ["XO", "GAMMA"]:
                        rdata_lat.append(fr"${val:.2f} \pm {unc:.2f}$")

                otab.add_row(rdata)
                otab_lat.add_row(rdata_lat)

        otab.write(
            f"tables/gor24_smc{ctype}_ensemble_params.dat", format="ascii.ipac", overwrite=True
        )

        otab_lat.write(
            f"tables/gor24_smc{ctype}_ensemble_dust_params.tex",
            format="aastex",
            col_align="lcccc",
            latexdict={
                "caption": r"Extinction Parameters \label{tab_ext_col_param}",
            },
            overwrite=True,
        )