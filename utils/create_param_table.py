import glob

# import numpy as np
from astropy.table import QTable

# import astropy.units as u

from measure_extinction.extdata import ExtData

# read in the Valencic et al. 2004 data and write the common format table
if __name__ == "__main__":

    for ctype in ["", "_forecor"]:
        # get all the files with FM paramters
        files = glob.glob(f"fits/*ext{ctype}_FM90.fits")
        n_files = len(files)

        otab = QTable(
            names=("name", "AV", "EBV", "RV", "C1", "C2", "C3", "C4", "x0", "gamma"),
            dtype=("S", "f", "f", "f", "f", "f", "f", "f", "f", "f"),
        )
        for cfile in files:
            edata = ExtData(filename=cfile)
            edata.calc_RV()

            rdata = []
            rdata.append(cfile)
            rdata.append(edata.columns["AV"][0])
            rdata.append(edata.columns["EBV"][0])
            rv = edata.columns["RV"][0]
            rdata.append(rv)
            rdata.append((edata.fm90_p50_fit["C1"][0] - 1.0) * rv)
            rdata.append(edata.fm90_p50_fit["C2"][0] * rv)
            rdata.append(edata.fm90_p50_fit["C3"][0] * rv)
            rdata.append(edata.fm90_p50_fit["C4"][0] * rv)
            rdata.append(edata.fm90_p50_fit["XO"][0])
            rdata.append(edata.fm90_p50_fit["GAMMA"][0])

            otab.add_row(rdata)

        otab.write(f"gor24_smc{ctype}_ensemble_params.dat", format="ascii.ipac", overwrite=True)
