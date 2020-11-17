import matplotlib.pyplot as plt
import argparse
import matplotlib

from plot_extinction import plot_ext_stack

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--eps", help="save figure as an eps file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    fontsize = 14

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    fig, ax = plt.subplots(figsize=(8, 6))

    datapath = "/home/kgordon/Hubble/SMCExt/SENDTOKARL_FINAL/"
    plot_ext_stack(
        "data/smc_stars_jwstprop.dat",
        ax,
        idlsave=True,
        locpath=datapath,
        fontsize=fontsize,
    )
    # ax[0].set_xlim(0.0, 9.0)

    ax.set_xlabel(r"$1/\lambda$ [$\mu m^{-1}$]", fontsize=1.3 * fontsize)
    ax.set_ylabel(r"$E(\lambda - V)/E(B - V)$ + offset", fontsize=1.3 * fontsize)

    fig.tight_layout()

    if args.png:
        fig.savefig("smc_ext_jwstprop.png")
    elif args.pdf:
        fig.savefig("smc_ext_jwstprop.pdf")
    elif args.eps:
        fig.savefig("smc_ext_jwstprop.eps")
    else:
        plt.show()
