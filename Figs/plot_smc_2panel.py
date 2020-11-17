import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
import argparse
import matplotlib

from plot_extinction import plot_ext_stack

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--eps", help="save figure as an eps file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    fontsize = 16

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    fig, ax = pyplot.subplots(figsize=(15, 9))
    gs = gridspec.GridSpec(2, 2)
    ax = []
    ax.append(pyplot.subplot(gs[0:2, 0]))
    ax.append(pyplot.subplot(gs[0:2, 1]))
    # ax.append(pyplot.subplot(gs[1, 1]))

#    datapath = "/home/kgordon/Hubble/SMCExt/Ed/"
    datapath = "/home/kgordon/Hubble/SMCExt/SENDTOKARL_FINAL/"
    plot_ext_stack(
        "data/smc_stars_reddened_good_highebv.dat",
        ax[0],
        idlsave=True,
        locpath=datapath,
        fontsize=fontsize,
    )
    ax[0].set_xlim(0.0, 9.0)
    ylimits = ax[0].get_ylim()
    xlimits = ax[0].get_xlim()
    ax[0].text(
        xlimits[0] + 0.05 * (xlimits[1] - xlimits[0]),
        ylimits[0] + 0.95 * (ylimits[1] - ylimits[0]),
        u"SMC, Good, $E(B-V) > 0.15$",
        fontsize=1.5 * fontsize,
        horizontalalignment="left",
    )

    plot_ext_stack(
        "data/smc_stars_reddened_good_lowebv.dat",
        ax[1],
        idlsave=True,
        locpath=datapath,
        fontsize=fontsize,
    )
    ax[1].set_xlim(0.0, 9.0)
    ylimits = ax[1].get_ylim()
    xlimits = ax[1].get_xlim()
    ax[1].text(
        xlimits[0] + 0.05 * (xlimits[1] - xlimits[0]),
        ylimits[0] + 0.9 * (ylimits[1] - ylimits[0]),
        u"Foreground, Good, $E(B-V) < 0.15$",
        fontsize=fontsize,
        horizontalalignment="left",
    )

    ax[0].set_xlabel(r"$1/\lambda$ [$\mu m^{-1}$]", fontsize=1.3 * fontsize)
    ax[0].set_ylabel(r"$E(\lambda - V)/E(B - V)$ + offset", fontsize=1.3 * fontsize)

    ax[0].set_ylim(-5, 100.)
    ax[1].set_ylim(-5, 100.)

    # ax[2].set_xlabel(r"$1/\lambda$ [$\mu m^{-1}$]", fontsize=1.3 * fontsize)

    fig.tight_layout()

    if args.png:
        fig.savefig("smc_ext_2panel.png")
    elif args.pdf:
        fig.savefig("smc_ext_2panel.pdf")
    elif args.eps:
        fig.savefig("smc_ext_2panel.eps")
    else:
        pyplot.show()
