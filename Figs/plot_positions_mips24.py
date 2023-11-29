import argparse
import matplotlib.pyplot as plt

from astropy.wcs import WCS
from astropy.io import fits
from astropy.table import QTable
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.visualization import SqrtStretch, ImageNormalize


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    filename = '~/Spitzer/SAGE-SMC/delivery/smc_e0_e1_e2_24_all_20jan09.cal.delivery.fits'

    hdu = fits.open(filename)[0]
    data = hdu.data
    wcs = WCS(hdu.header)

    # make a smaller image
    coord = SkyCoord("01:00:00", "-72:45:00", unit=(u.hourangle, u.deg))
    cutout = Cutout2D(data, coord, (3.0 * u.deg, 4.0 * u.deg), wcs=wcs)
    data = cutout.data
    wcs = cutout.wcs

    ax = plt.subplot(projection=wcs)
    fig = plt.gcf()
    fig.set_size_inches(8.5, 6)

    norm = ImageNormalize(vmin=0.0, vmax=5.0, stretch=SqrtStretch())

    # ax.set_figsize((10, 6))
    ax.imshow(data, vmin=-0.1, vmax=2, origin='lower', cmap="binary")
    ax.imshow(data, origin='lower', norm=norm, cmap="binary")
    ax.grid(color='black', ls='dotted', alpha=0.5)
    ax.set_xlabel("DEC")
    ax.set_ylabel("RA")

    # read in the nobump stars and overplot
    ptab = QTable.read("data/smc_positions_nobump.dat", format="ascii.commented_header")
    for k in range(len(ptab)):
        coord = SkyCoord(
            ptab["ra"][k],
            ptab["dec"][k],
            unit=(u.hourangle, u.deg),
        )
        ax.scatter(coord.ra.degree, coord.dec.degree, transform=ax.get_transform('fk5'), s=75,
                   edgecolor='blue', facecolor='none', linewidth=2, alpha=0.75)
        ax.annotate(ptab["name"][k], wcs.world_to_pixel(coord), color="b")

    # read in the bump stars and overplot
    ptab = QTable.read("data/smc_positions_bump.dat", format="ascii.commented_header")
    for k in range(len(ptab)):
        coord = SkyCoord(
            ptab["ra"][k],
            ptab["dec"][k],
            unit=(u.hourangle, u.deg),
        )
        ax.scatter(coord.ra.degree, coord.dec.degree, transform=ax.get_transform('fk5'), s=75,
                   edgecolor='red', facecolor='none', linewidth=2, alpha=0.75, marker="P")
        ax.annotate(ptab["name"][k], wcs.world_to_pixel(coord), color="r")

    #ax.scatter(266.78238, -28.769255, transform=ax.get_transform('fk5'), s=300,
    #           edgecolor='white', facecolor='none')

    plt.tight_layout()

    save_str = "smc_mips24_positions"
    if args.png:
        plt.savefig(f"{save_str}.png")
    elif args.pdf:
        plt.savefig(f"{save_str}.pdf")
    else:
        plt.show()