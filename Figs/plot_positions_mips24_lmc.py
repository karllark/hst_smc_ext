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

    #filename = '~/Spitzer/SAGE-SMC/delivery/smc_e0_e1_e2_24_all_20jan09.cal.delivery.fits'
    filename_qpah = "data/lmc_dl07_fixedbeta_06sep18_qPAH.fits.gz"
    filename = filename_qpah

    hdu = fits.open(filename)[0]
    data = hdu.data
    wcs = WCS(hdu.header)

    # make a smaller image
    coord = SkyCoord("05:24:00", "-68:45:00", unit=(u.hourangle, u.deg))
    cutout = Cutout2D(data, coord, (8.0 * u.deg, 8.0 * u.deg), wcs=wcs)
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
    fig.set_size_inches(8.5, 8.25)

    norm = ImageNormalize(vmin=0.0, vmax=3.0, stretch=SqrtStretch())

    # ax.set_figsize((10, 6))
    ax.imshow(data, vmin=-0.1, vmax=2, origin='lower', cmap="binary")
    ax.imshow(data, origin='lower', norm=norm)  #, cmap="binary")
    ax.grid(color='black', ls='dotted', alpha=0.5)
    ax.set_xlabel("RA")
    ax.set_ylabel("DEC")

    fnames = ["avg", "lmc2"]
    fsyms = ["o", "s"]
    fcols = ["b", "r"]

    names = []
    qpahs = []
    qpahs_unc = []
    for cname, csym, ccol in zip(fnames, fsyms, fcols):
        ptab = QTable.read(f"data/lmc_positions_{cname}.dat", format="ascii.commented_header")
        for k in range(len(ptab)):
            coord = SkyCoord(
                ptab["ra"][k],
                ptab["dec"][k],
                unit=(u.hourangle, u.deg),
            )
            ax.scatter(coord.ra.degree, coord.dec.degree, transform=ax.get_transform('fk5'), s=45,
                    edgecolor=ccol, facecolor='none', linewidth=2, alpha=0.75, marker=csym)
            if args.names:
                ax.annotate(ptab["name"][k], wcs.world_to_pixel(coord), color=ccol)
            # get the qpah values from the map
            pix_qpah = wcs_qpah.world_to_pixel(coord)
            if ((pix_qpah[1] >= 0)
                & (pix_qpah[1] < data_qpah[0].shape[0] - 1)
                & (pix_qpah[0] >= 0)
                & (pix_qpah[0] < data_qpah[0].shape[1] - 1)):
                names.append(ptab["name"].data[k])
                qpahs.append(data_qpah[0][int(pix_qpah[1]), int(pix_qpah[0])])
                qpahs_unc.append(0.5 * (data_qpah[1][int(pix_qpah[1]), int(pix_qpah[0])] - data_qpah[2][int(pix_qpah[1]), int(pix_qpah[0])]))
            else:
                print(ptab["name"].data[k], " not in qpah image")

    # output qpah values
    otab = QTable()
    otab["name"] = names
    otab["qpah"] = qpahs
    otab["qpah_unc"] = qpahs_unc
    otab.write("tables/qpah_from_image_lmc.dat", format="ascii.ipac", overwrite=True)

    legend_elements = [Line2D([0], [0], marker='o', markerfacecolor="none", color="none", label='LMC Average Sample', 
                              markeredgecolor='b', markersize=10),
                       Line2D([0], [0], marker='s', markerfacecolor="none", color="none", label='LMC2 Supershell Sample', 
                              markeredgecolor='r', markersize=10)]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=0.8*fontsize)

    ax.annotate(r"QPAH", (100, 100))

    plt.tight_layout()

    save_str = "lmc_qpah_positions"
    if args.png:
        plt.savefig(f"{save_str}.png")
    elif args.pdf:
        plt.savefig(f"{save_str}.pdf")
    else:
        plt.show()