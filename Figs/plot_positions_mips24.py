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

    filename = '~/Spitzer/SAGE-SMC/delivery/smc_e0_e1_e2_24_all_20jan09.cal.delivery.fits'

    hdu = fits.open(filename)[0]
    data = hdu.data
    wcs = WCS(hdu.header)

    # make a smaller image
    coord = SkyCoord("01:11:00", "-72:40:00", unit=(u.hourangle, u.deg))
    cutout = Cutout2D(data, coord, (3.0 * u.deg, 5.0 * u.deg), wcs=wcs)
    data = cutout.data
    wcs = cutout.wcs

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
    fig.set_size_inches(8.5, 6)

    norm = ImageNormalize(vmin=0.0, vmax=3.0, stretch=SqrtStretch())

    # ax.set_figsize((10, 6))
    ax.imshow(data, vmin=-0.1, vmax=2, origin='lower', cmap="binary")
    ax.imshow(data, origin='lower', norm=norm, cmap="binary")
    ax.grid(color='black', ls='dotted', alpha=0.5)
    ax.set_xlabel("RA")
    ax.set_ylabel("DEC")

    # read in the low ebv stars and overplot
    ptab = QTable.read("data/smc_positions_lowebv.dat", format="ascii.commented_header")
    for k in range(len(ptab)):
        coord = SkyCoord(
            ptab["ra"][k],
            ptab["dec"][k],
            unit=(u.hourangle, u.deg),
        )
        ax.scatter(coord.ra.degree, coord.dec.degree, transform=ax.get_transform('fk5'), s=45,
                   edgecolor='tab:brown', facecolor='none', linewidth=2, alpha=0.75, marker="v")
        if args.names:
            ax.annotate(ptab["name"][k], wcs.world_to_pixel(coord), color="tab:brown")

    # read in the nobump stars and overplot
    ptab = QTable.read("data/smc_positions_nobump.dat", format="ascii.commented_header")
    for k in range(len(ptab)):
        coord = SkyCoord(
            ptab["ra"][k],
            ptab["dec"][k],
            unit=(u.hourangle, u.deg),
        )
        ax.scatter(coord.ra.degree, coord.dec.degree, transform=ax.get_transform('fk5'), s=65,
                   edgecolor='blue', facecolor='none', linewidth=2, alpha=0.75)
        if args.names:
            ax.annotate(ptab["name"][k], wcs.world_to_pixel(coord), color="b")

    # read in the flat stars and overplot
    ptab = QTable.read("data/smc_positions_flat.dat", format="ascii.commented_header")
    for k in range(len(ptab)):
        coord = SkyCoord(
            ptab["ra"][k],
            ptab["dec"][k],
            unit=(u.hourangle, u.deg),
        )
        ax.scatter(coord.ra.degree, coord.dec.degree, transform=ax.get_transform('fk5'), s=55,
                   edgecolor='cyan', facecolor='none', linewidth=2, alpha=0.75, marker="s")
        if args.names:
            ax.annotate(ptab["name"][k], wcs.world_to_pixel(coord), color="c")

    # read in the bump stars and overplot
    ptab = QTable.read("data/smc_positions_bump.dat", format="ascii.commented_header")
    for k in range(len(ptab)):
        coord = SkyCoord(
            ptab["ra"][k],
            ptab["dec"][k],
            unit=(u.hourangle, u.deg),
        )
        ax.scatter(coord.ra.degree, coord.dec.degree, transform=ax.get_transform('fk5'), s=65,
                   edgecolor='red', facecolor='none', linewidth=2, alpha=0.75, marker="P")
        if args.names:
            ax.annotate(ptab["name"][k], wcs.world_to_pixel(coord), color="r")

    legend_elements = [Line2D([0], [0], marker='o', markerfacecolor="none", color="none", label='Weak/absent 2175 A bump', 
                              markeredgecolor='b', markersize=10),
                       Line2D([0], [0], marker='P', markerfacecolor="none", color="none", label='Significant 2175 A bump', 
                              markeredgecolor='r', markersize=10),
                       Line2D([0], [0], marker='s', markerfacecolor="none", color="none", label='Flat', 
                              markeredgecolor='c', markersize=10),
                       Line2D([0], [0], marker='v', markerfacecolor="none", color="none", label=r'$E(B-V)_\mathrm{SMC} < 0.1$', 
                              markeredgecolor='tab:brown', markersize=10)]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=0.8*fontsize)

    ax.annotate(r"MIPS 24 $\mu$m", (100, 100))

    plt.tight_layout()

    save_str = "smc_mips24_positions"
    if args.png:
        plt.savefig(f"{save_str}.png")
    elif args.pdf:
        plt.savefig(f"{save_str}.pdf")
    else:
        plt.show()