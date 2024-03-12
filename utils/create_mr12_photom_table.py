from astropy.table import QTable


if __name__ == "__main__":
    # read in the SMIDGE BEAST catalog for the MR12 stars
    # created by Petia Yanacholova Jones
    tab = QTable.read("MR12_ext/SMIDGE_BEAST_MR12_sorted.fits")

    table_lat = QTable(
        names=(
            "Name",
            "F336W", "F475W", "F550M", "F814W", "F110W", "F160W",
        ),
        dtype=("S", "S", "S", "S", "S", "S", "S"),
    )

    names = ["MR12-09", "MR12-10", "MR12-11"]
    bands = ["F336W", "F475W", "F550M", "F814W", "F110W", "F160W"]

    for k, cname in enumerate(names):
        rdata_lat = []
        rdata_lat.append(cname)
        rownum = k + 8
        for cband in bands:
            val = tab[f"{cband}_VEGA"][rownum]
            unc = tab[f"{cband}_ERR"][rownum]
            rdata_lat.append(f"${val:.3f} \pm {unc:.3f}$")

        table_lat.add_row(rdata_lat)

    table_lat.write(
        "tables/mr12_phot.tex",
        format="aastex",
        col_align="lcccccc",
        latexdict={
            "caption": r"Sample HST photometry \label{tab_hst}",
        },
        overwrite=True,
    )