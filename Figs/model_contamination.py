import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import astropy.units as u

from dust_extinction.parameter_averages import G16


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # plots
    fontsize = 14

    font = {"size": fontsize}

    plt.rc("font", **font)

    plt.rc("lines", linewidth=1)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("xtick.minor", width=2)
    plt.rc("ytick.major", width=2)
    plt.rc("ytick.minor", width=2)

    figsize = (14.0, 5.0)
    fig, ax = plt.subplots(ncols=2, figsize=figsize)

    mod_x = np.logspace(np.log10(0.11), np.log10(2.5), num=1000) * u.micron
    ebvtots = [0.15, 0.30]
    ebvfore1 = 0.04
    ebvfore2 = 0.08

    RvA = 3.1
    G16_lowebv1 = G16(RvA=RvA, fA=1.0)
    G16_lowebv0 = G16(RvA=RvA, fA=0.0)

    for k, ebvtot in enumerate(ebvtots):
        G16_lowebv_fore1 = G16(RvA=RvA, fA=ebvfore1 / ebvtot)
        G16_lowebv_fore2 = G16(RvA=RvA, fA=ebvfore2 / ebvtot)
        ax[k].plot(
            1.0 / mod_x,
            G16_lowebv0(mod_x),
            "k-.",
            lw=2,
            alpha=0.5,
            label=r"$E(B-V)_{MW} = 0.00$",
        )
        ax[k].plot(
            1.0 / mod_x,
            G16_lowebv_fore1(mod_x),
            "b-",
            lw=2,
            alpha=0.75,
            label=r"$E(B-V)_{MW} = 0.04$",
        )
        ax[k].plot(
            1.0 / mod_x,
            G16_lowebv_fore2(mod_x),
            "g-",
            lw=2,
            alpha=0.75,
            label=r"$E(B-V)_{MW} = 0.08$",
        )
        ax[k].plot(
            1.0 / mod_x,
            G16_lowebv1(mod_x),
            "k:",
            lw=2,
            alpha=0.5,
            label=r"$E(B-V)_{MW} = " + f"{ebvtot}" + "$",
        )

        ax[k].set_xlim(0.0, 9.0)
        ax[k].set_xlabel(r"$1/\lambda$ $[\mu m^{-1}]$")
        ax[k].set_ylabel(r"$A(\lambda)/A(V)$")
        ax[k].legend(title=r"$E(B-V)_{total} = " + f"{ebvtot}" + "$")

    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")
    ax[1].xaxis.set_major_formatter(ScalarFormatter())
    ax[1].yaxis.set_major_formatter(ScalarFormatter())

    fig.tight_layout()

    sfile = "fore_cont_model"
    if args.png:
        fig.savefig(f"{sfile}.png")
    elif args.pdf:
        fig.savefig(f"{sfile}.pdf")
    else:
        plt.show()
