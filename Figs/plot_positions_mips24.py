import argparse
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from astropy.wcs import WCS
from astropy.io import fits
from astropy.table import QTable
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.visualization import SqrtStretch, ImageNormalize


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--names", help="annotate with star names", action="store_true")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    filename = (
        "~/Spitzer/SAGE-SMC/delivery/smc_e0_e1_e2_24_all_20jan09.cal.delivery.fits"
    )
    filename_qpah = "data/smc_dl07_fixedbeta_06sep18_1sig_qPAH.fits.gz"
    # filename = filename_qpah

    hdu = fits.open(filename)[0]
    data = hdu.data
    wcs = WCS(hdu.header)

    # make a smaller image
    coord = SkyCoord("01:11:00", "-72:40:00", unit=(u.hourangle, u.deg))
    cutout = Cutout2D(data, coord, (3.0 * u.deg, 5.0 * u.deg), wcs=wcs)
    data = cutout.data
    wcs = cutout.wcs

    # qpah images
    data_qpah = []
    hdu_qpah = fits.open(filename_qpah)[0]
    data_qpah.append(hdu_qpah.data)
    wcs_qpah = WCS(hdu_qpah.header)
    # p16 and p84 for uncs
    hdu_qpah_p = fits.open(filename_qpah.replace("qPAH", "qPAH_84th"))[0]
    data_qpah.append(hdu_qpah_p.data)
    hdu_qpah_m = fits.open(filename_qpah.replace("qPAH", "qPAH_16th"))[0]
    data_qpah.append(hdu_qpah_m.data)

    fontsize = 14

    font = {"size": fontsize}

    plt.rc("font", **font)

    plt.rc("lines", linewidth=1)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("xtick.minor", width=2)
    plt.rc("ytick.major", width=2)
    plt.rc("ytick.minor", width=2)

    ax = plt.subplot(projection=wcs)
    fig = plt.gcf()
    fig.set_size_inches(8.5, 5.25)

    norm = ImageNormalize(vmin=0.0, vmax=3.0, stretch=SqrtStretch())
    norm2 = ImageNormalize(vmin=0.0, vmax=20.0, stretch=SqrtStretch())

    # ax.set_figsize((10, 6))
    # ax.imshow(data, vmin=-0.1, vmax=2, origin='lower', cmap="binary")
    ax.imshow(data, origin="lower", norm=norm, cmap="binary")
    ax.grid(color="black", ls="dotted", alpha=0.5)
    ax.set_xlabel("RA")
    ax.set_ylabel("DEC")

    # add zoom region
    pixscale = 0.000691667 * 3600.0  # in arcsec/pixel - from fits header
    extent = (-3, 4, -4, 3)
    cx, cy = 6145.0, 865.0
    chwidth = 10
    x1, x2, y1, y2 = (
        cx - chwidth,
        cx + chwidth,
        cy - chwidth,
        cy + chwidth,
    )  # subregion of the original image
    axins = ax.inset_axes(
        [0.8, 0.7, 0.2, 0.2],
        xlim=(x1, x2),
        ylim=(y1, y2),
        xticklabels=[],
        yticklabels=[],
    )
    # axins.imshow(data, extent=extent, origin="lower")
    axins.imshow(data, origin="lower", norm=norm2, cmap="binary")
    ax.indicate_inset_zoom(axins, edgecolor="black")
    # add scale
    ax.annotate(
        f'{2.*pixscale*chwidth:.0f}"', (6400.0, 4000.0), fontsize=0.8 * fontsize
    )

    # add 2nd zoom region
    extent = (-3, 4, -4, 3)
    cx, cy = 3165.0, 1315.0
    chwidth = 40.0
    x1, x2, y1, y2 = (
        cx - chwidth,
        cx + chwidth,
        cy - chwidth,
        cy + chwidth,
    )  # subregion of the original image
    axins2 = ax.inset_axes(
        [0.1, 0.45, 0.2, 0.2],
        xlim=(x1, x2),
        ylim=(y1, y2),
        xticklabels=[],
        yticklabels=[],
    )
    axins2.imshow(data, origin="lower", norm=norm2, cmap="binary")
    ax.indicate_inset_zoom(axins2, edgecolor="black")
    # add scale
    ax.annotate(
        f'{2.*pixscale*chwidth:.0f}"', (1300.0, 2900.0), fontsize=0.8 * fontsize
    )

    fnames = ["lowebv", "nobump", "flat", "bump"]
    fsyms = ["v", "o", "s", "P"]
    fcols = ["tab:brown", "b", "c", "r"]

    names = []
    qpahs = []
    qpahs_unc = []
    for cname, csym, ccol in zip(fnames, fsyms, fcols):
        ptab = QTable.read(
            f"data/smc_positions_{cname}.dat", format="ascii.commented_header"
        )
        for k in range(len(ptab)):
            coord = SkyCoord(
                ptab["ra"][k],
                ptab["dec"][k],
                unit=(u.hourangle, u.deg),
            )
            ax.scatter(
                coord.ra.degree,
                coord.dec.degree,
                transform=ax.get_transform("fk5"),
                s=45,
                edgecolor=ccol,
                facecolor="none",
                linewidth=2,
                alpha=0.75,
                marker=csym,
            )

            # zooms
            cploc = wcs.world_to_pixel(coord)
            axins.scatter(
                [cploc[0]],
                [cploc[1]],
                s=45,
                edgecolor=ccol,
                facecolor="none",
                linewidth=2,
                alpha=0.75,
                marker=csym,
            )
            axins2.scatter(
                [cploc[0]],
                [cploc[1]],
                s=45,
                edgecolor=ccol,
                facecolor="none",
                linewidth=2,
                alpha=0.75,
                marker=csym,
            )

            # print(ptab["name"][k], wcs.world_to_pixel(coord))
            if args.names:
                ax.annotate(ptab["name"][k], wcs.world_to_pixel(coord), color=ccol)
            # get the qpah values from the map
            pix_qpah = wcs_qpah.world_to_pixel(coord)
            if (
                (pix_qpah[1] >= 0)
                & (pix_qpah[1] < data_qpah[0].shape[0] - 1)
                & (pix_qpah[0] >= 0)
                & (pix_qpah[0] < data_qpah[0].shape[1] - 1)
            ):
                names.append(ptab["name"].data[k])
                qpahs.append(data_qpah[0][int(pix_qpah[1]), int(pix_qpah[0])])
                qpahs_unc.append(
                    0.5
                    * (
                        data_qpah[1][int(pix_qpah[1]), int(pix_qpah[0])]
                        - data_qpah[2][int(pix_qpah[1]), int(pix_qpah[0])]
                    )
                )
            else:
                print(ptab["name"].data[k], " not in qpah image")

    # output qpah values
    otab = QTable()
    otab["name"] = names
    otab["qpah"] = qpahs
    otab["qpah_unc"] = qpahs_unc
    otab.write("tables/qpah_from_image.dat", format="ascii.ipac", overwrite=True)

    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            markerfacecolor="none",
            color="none",
            label="Weak/absent 2175 A bump",
            markeredgecolor="b",
            markersize=10,
        ),
        Line2D(
            [0],
            [0],
            marker="P",
            markerfacecolor="none",
            color="none",
            label="Significant 2175 A bump",
            markeredgecolor="r",
            markersize=10,
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            markerfacecolor="none",
            color="none",
            label="Flat",
            markeredgecolor="c",
            markersize=10,
        ),
        Line2D(
            [0],
            [0],
            marker="v",
            markerfacecolor="none",
            color="none",
            label=r"$E(B-V)_\mathrm{SMC} < 0.1$",
            markeredgecolor="tab:brown",
            markersize=10,
        ),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=0.8 * fontsize)

    ax.annotate(r"MIPS 24 $\mu$m", (100, 100))

    plt.tight_layout()

    save_str = "smc_mips24_positions"
    if args.png:
        plt.savefig(f"{save_str}.png")
    elif args.pdf:
        plt.savefig(f"{save_str}.pdf")
    else:
        plt.show()
