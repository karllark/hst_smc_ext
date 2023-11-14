import glob

import numpy as np
from astropy.table import QTable

# import astropy.units as u

from measure_extinction.extdata import ExtData

# read in the Valencic et al. 2004 data and write the common format table
if __name__ == "__main__":

    for ctype in ["_forecor"]:
        file1 = "data/smc_stars_reddened_good.dat"

        f = open(file1, "r")
        file_lines = list(f)
        starnames = []
        for line in file_lines:
            if (line.find("#") != 0) & (len(line) > 0):
                name = line.rstrip()
                starnames.append(name)
        starnames = np.sort(starnames)

        otab = QTable(
            names=("name", "AV", "AV_unc", "EBV", "EBV_unc", "RV", "RV_unc",
                   "C1", "C1_unc", "C2", "C2_unc", "B3", "B3_unc", "C4", "C4_unc",
                   "x0", "x0_unc", "gamma", "gamma_unc"),
            dtype=("S", "f", "f", "f", "f", "f", "f", "f", "f", "f", "f", "f",
                   "f", "f", "f", "f", "f", "f", "f"),
        )
        for cname in starnames:
            cfile = f"fits/{cname}_ext_forecor_FM90.fits"

            edata = ExtData(filename=cfile)
            edata.calc_RV()

            rdata = []
            rdata.append(cname)
            rdata.append(edata.columns["AV"][0])
            rdata.append(edata.columns["AV"][1])
            rdata.append(edata.columns["EBV"][0])
            rdata.append(edata.columns["EBV"][1])
            rdata.append(edata.columns["RV"][0])
            rdata.append(edata.columns["RV"][1])
            rdata.append(edata.fm90_p50_fit["C1"][0])
            rdata.append(edata.fm90_p50_fit["C1"][1])
            rdata.append(edata.fm90_p50_fit["C2"][0])
            rdata.append(edata.fm90_p50_fit["C2"][1])
            rdata.append(edata.fm90_p50_fit["B3"][0])
            rdata.append(edata.fm90_p50_fit["B3"][1])
            rdata.append(edata.fm90_p50_fit["C4"][0])
            rdata.append(edata.fm90_p50_fit["C4"][1])
            rdata.append(edata.fm90_p50_fit["XO"][0])
            rdata.append(edata.fm90_p50_fit["XO"][1])
            rdata.append(edata.fm90_p50_fit["GAMMA"][0])
            rdata.append(edata.fm90_p50_fit["GAMMA"][1])

            otab.add_row(rdata)

        otab.write(f"gor24_smc{ctype}_ensemble_params.dat", format="ascii.ipac", overwrite=True)
