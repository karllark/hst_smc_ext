import glob

import numpy as np
from astropy.table import QTable

# import astropy.units as u

from helpers import prettyname


if __name__ == "__main__":

    # read in the spreadsheet data
    sdata = QTable.read("data/edf_fits_param_edit.csv")
    print(sdata.colnames)

    filename = "data/smc_stars_all.dat"
    f = open(filename, "r")
    file_lines = list(f)
    extnames = []
    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            extnames.append(name)

    pnames = []
    for k, cname in enumerate(extnames):
        pname = prettyname(cname)
        pnames.append(pname)

        print(pname.replace(" ", "_"))
        (mindx,) = np.where(pname.replace(" ", "_") == sdata["name"])
        if len(mindx) == 0:
            print(f"{cname} not found in params file")