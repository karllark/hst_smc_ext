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

    table_lat = QTable(
        names=(
            "Name",
            r"$T_\mathrm{eff}$ [K]",
            r"$\log (g)$",
        ),
        dtype=("str", "str", "str"),
    )

    pnames = []
    for k, cname in enumerate(extnames):
        pname = prettyname(cname)
        pnames.append(pname)

        print(pname.replace(" ", "_"))
        (mindx,) = np.where(pname.replace(" ", "_") == sdata["name"])
        if len(mindx) == 0:
            print(f"{cname} not found in params file")
        else:
            teff = sdata["log T"][mindx].data[0]
            teff_unc = np.sqrt(
                sdata["sig"][mindx].data[0] ** 2 + sdata["sig(MC)"][mindx].data[0] ** 2
            )
            pteff = 10 ** (teff + teff_unc)
            mteff = 10 ** (teff - teff_unc)
            teff = 0.5 * (pteff + mteff)
            teff_unc = 0.5 * (pteff - mteff)
            logg = sdata["logg"][mindx].data[0]
            logg_unc = np.sqrt(
                sdata["logg_sig"][mindx].data[0] ** 2
                + sdata["logg_sigmc"][mindx].data[0] ** 2
            )
            table_lat.add_row(
                (
                    pname,
                    f"{teff:.0f} +/- {teff_unc:.0f}",
                    f"{logg:.3f} +/- {logg_unc:.3f}",
                )
            )

    table_lat.write(
        "tables/stellar_params.tex",
        format="aastex",
        col_align="lccc",
        latexdict={
            "caption": r"Stellar Fit Parameters \label{tab_stellar_param}",
        },
        overwrite=True,
    )
