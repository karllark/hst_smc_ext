from astropy.table import QTable


if __name__ == "__main__":
    # read in the SMIDGE BEAST catalog for the MR12 stars
    # created by Petia Yanacholova Jones
    tab = QTable.read("MR12_ext/SMIDGE_BEAST_MR12_sorted.fits")

    # table rotated from others as too many bands
    table_lat = QTable(
        names=(
            "Band",
            "MR12 09", "MR12 10", "MR12 11",
        ),
        dtype=("S", "S", "S", "S"),
    )

    names = ["MR12-09", "MR12-10", "MR12-11"]
    bands = ["F225W", "F275W", "F336W", "F475W", "F550M", "F814W", "F110W", "F160W"]

    for cband in bands:
        rdata_lat = []
        rdata_lat.append(cband)
        for k, cname in enumerate(names):
            rownum = k + 8
            val = tab[f"{cband}_VEGA"][rownum]
            unc = tab[f"{cband}_ERR"][rownum]
            rdata_lat.append(f"${val:.3f} \pm {unc:.3f}$")

        table_lat.add_row(rdata_lat)

    table_lat.write(
        "tables/mr12_phot.tex",
        format="aastex",
        col_align="lccc",
        latexdict={
            "caption": r"Sample HST photometry \label{tab_hst}",
        },
        overwrite=True,
    )