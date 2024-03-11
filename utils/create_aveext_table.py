import numpy as np
from astropy.table import QTable

if __name__ == "__main__":

    otab_lat = QTable(
        # fmt: off
        names=("$\lambda$", "$x$", "Band", "Average", "Bumps"),
        dtype=("S", "S", "S", "S", "S")
        # fmt:on
    )

    avetab = QTable.read("SMC_Average_Gordon24_ext.dat", format="ascii.commented_header")
    bumptab = QTable.read("SMC_Bumps_Gordon24_ext.dat", format="ascii.commented_header")

    sindxs, = np.where(avetab["type"] == "spec")
    for k in sindxs:
        rdata_lat = []
        rdata_lat.append(fr"${avetab['wave'][k]:.3f}$")
        rdata_lat.append(fr"${1./avetab['wave'][k]:.3f}$")
        rdata_lat.append("IUE/STIS")
        rdata_lat.append(fr"${avetab['A(l)/A(V)'][k]:.3f} \pm {avetab['unc'][k]:.3f}$")
        rdata_lat.append(fr"${bumptab['A(l)/A(V)'][k]:.3f} \pm {bumptab['unc'][k]:.3f}$")
        otab_lat.add_row(rdata_lat)

    #bnames = ["WFC3_F336W", "JohnU", "JohnB", "ACS_F475W", "JohnV", "ACS_F550M",
    #          "ACS_F814W", "WFC3_F110W", "JohnJ", "WFC3_F160W", "JohnH", "JohnK"]
    bnames = ["JohnU", "JohnB", "JohnV", "JohnJ", "JohnH", "JohnK"]
    for cname in bnames:
        rdata_lat = []
        k2, = np.where(avetab["type"] == cname)
        k1, = np.where(bumptab["type"] == cname)
        if len(k2) > 0:
            k2 = k2[0]
        else:
            k2 = None
        k1 = k1[0]

        rdata_lat.append(fr"${bumptab['wave'][k1]:.3f}$")
        rdata_lat.append(fr"${1./bumptab['wave'][k1]:.3f}$")
        rdata_lat.append(cname.replace("_", " "))
        if k2 is not None:
            rdata_lat.append(fr"${avetab['A(l)/A(V)'][k2]:.3f} \pm {avetab['unc'][k2]:.3f}$")
        else:
            rdata_lat.append("\\nodata")
        rdata_lat.append(fr"${bumptab['A(l)/A(V)'][k1]:.3f} \pm {bumptab['unc'][k1]:.3f}$")

        otab_lat.add_row(rdata_lat)

    otab_lat.write(
        f"tables/gor24_average_curves.tex",
        format="aastex",
        col_align="lcccc",
        latexdict={
            "caption": r"Sample Averages \label{tab_aveext}",
        },
        overwrite=True,
    )