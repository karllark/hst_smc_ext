import matplotlib.pyplot as pyplot
# import matplotlib.gridspec as gridspec
import argparse
import matplotlib

from plot_extinction import plot_ext_stack

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--suspect", help="show suspect sightlines", action="store_true")
    parser.add_argument("--both", help="plot both good and adjusted samples", action="store_true")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    if args.both:
        file0 = "data/smc_stars_reddened_good_highebv_bumps.dat"
        file2 = "data/smc_stars_reddened_good_lowebv.dat"
        file1 = "data/smc_stars_reddened_good_highebv.dat"
        file3 = "data/smc_stars_reddened_good_highebv_flat.dat"
        forecor1 = True
        forecor2 = True
        ptitle1 = r"$E(B-V)_\mathrm{SMC} \geq 0.1$"
        ptitle2 = r"$E(B-V)_\mathrm{SMC} < 0.1$"
        figsize = (13, 10)
        textyval = 0.90
        adjusted = False
    else:
        if args.suspect:
            file1 = "data/smc_stars_reddened_suspect.dat"
            # file3 = "data/smc_stars_reddened_suspect_adjusted.dat"
            stitle = "Set B"
            adjusted = False
            figsize = (13, 7)
            textyval = 0.90
        else:
            file1 = "data/smc_stars_reddened_good.dat"
            file1 = "data/smc_stars_all.dat"
            stitle = "Set A"
            adjusted = False
            figsize = (13, 13)
            textyval = 0.95
        file2 = file1
        forecor1 = False
        forecor2 = True
        ptitle1 = f"{stitle}, as measured"
        ptitle2 = f"{stitle}, foreground corrected"

    fontsize = 14

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    # fig, ax = pyplot.subplots(figsize=(13, 8))
    # gs = gridspec.GridSpec(2, 2)
    # ax = []
    # ax.append(pyplot.subplot(gs[0:2, 0]))
    # ax.append(pyplot.subplot(gs[0:2, 1]))
    # ax.append(pyplot.subplot(gs[1, 1]))
    #if args.suspect:
    #    fig, ax = pyplot.subplots(ncols=3, figsize=(19, 8))
    #else:
    fig, ax = pyplot.subplots(ncols=2, figsize=figsize, sharex=False, sharey=False)

    #    datapath = "/home/kgordon/Hubble/SMCExt/Ed/"
    datapath = "/home/kgordon/Python/hst_smc_ext/fits/"

    plot_ext_stack(
        file3,
        ax[0],
        locpath=datapath,
        fontsize=fontsize,
        forecor=forecor1,
        exoffset=0.
    )

    plot_ext_stack(
        file1,
        ax[0],
        locpath=datapath,
        fontsize=fontsize,
        forecor=forecor1,
        exoffset=15.
    )
    ax[0].set_xlim(0.0, 10.0)
    ax[0].set_ylim(-5, 110)
    ylimits = ax[0].get_ylim()
    xlimits = ax[0].get_xlim()
    # ax[0].text(
    #     xlimits[0] + 0.05 * (xlimits[1] - xlimits[0]),
    #     ylimits[0] + textyval * (ylimits[1] - ylimits[0]),
    #     ptitle1,
    #     fontsize=1.3 * fontsize,
    #     horizontalalignment="left",
    # )
    ax[0].text(1., 100., ptitle1, fontsize=1.1*fontsize, horizontalalignment="left")
    ax[0].text(9.25, 65.0, "Steep w/ Weak/Absent Bump", rotation=270., fontsize=0.9*fontsize,
               horizontalalignment="center", verticalalignment="center")
    ax[0].text(9.25, 10.0, "Flat", rotation=270., fontsize=0.9*fontsize,
               horizontalalignment="center", verticalalignment="center")

    plot_ext_stack(
        file0,
        ax[1],
        locpath=datapath,
        fontsize=fontsize,
        forecor=forecor1,
        exoffset=70.,
        topxaxis=False,
    )
    plot_ext_stack(
        file2,
        ax[1],
        locpath=datapath,
        fontsize=fontsize,
        forecor=forecor2,
        adjusted=adjusted,
    )
    ax[1].set_ylim(-5, 110)
    ax[1].set_xlim(0.0, 10.0)
    ylimits = ax[1].get_ylim()
    xlimits = ax[1].get_xlim()
    # ax[1].text(
    #     xlimits[0] + 0.05 * (xlimits[1] - xlimits[0]),
    #     ylimits[0] + textyval * (ylimits[1] - ylimits[0]),
    #     ptitle2,
    #     fontsize=1.3 * fontsize,
    #     horizontalalignment="left",
    # )
    ax[1].text(1., 100., ptitle1, fontsize=1.1*fontsize, horizontalalignment="left")
    ax[1].text(9.25, 90.0, "Significant Bump", rotation=270., fontsize=0.9*fontsize,
               horizontalalignment="center", verticalalignment="center")
    ax[1].text(1., 55., ptitle2, fontsize=1.1*fontsize, horizontalalignment="left")
    ax[1].text(9.25, 30.0, "Weak/Absent Bump", rotation=270., fontsize=0.9*fontsize,
               horizontalalignment="center", verticalalignment="center")

    # if args.suspect:
    #    plot_ext_stack(
    #        #  "data/smc_stars_reddened_good_lowebv.dat",
    #        file2,
    #        ax[2],
    #        locpath=datapath,
    #        fontsize=fontsize,
    #        forecor=forecor2,
    #        adjusted=True
    #    )
    #    ax[2].set_xlim(0.0, 9.0)
    #    ylimits = ax[2].get_ylim()
    #    xlimits = ax[2].get_xlim()
    #    ax[2].text(
    #        xlimits[0] + 0.05 * (xlimits[1] - xlimits[0]),
    #        ylimits[0] + 0.95 * (ylimits[1] - ylimits[0]),
    #        "Adjusted Foreground",
    #        fontsize=1.5 * fontsize,
    #        horizontalalignment="left",
    #    )

    ax[0].set_xlabel(r"$1/\lambda$ [$\mu m^{-1}$]", fontsize=1.3 * fontsize)
    ax[0].set_ylabel(r"$E(\lambda - V)/E(B - V)$ + offset", fontsize=1.3 * fontsize)

    # ax[0].set_ylim(-5, 100.0)
    # ax[1].set_ylim(-5, 100.0)

    ax[1].set_xlabel(r"$1/\lambda$ [$\mu m^{-1}$]", fontsize=1.3 * fontsize)
    ax[1].set_ylabel("")

    fig.tight_layout()

    filebase = "smc_mext"
    if args.suspect:
        filebase = f"{filebase}_suspect"
    if args.both:
        filebase = f"{filebase}_allext"
    if args.png:
        fig.savefig(f"{filebase}.png")
    elif args.pdf:
        fig.savefig(f"{filebase}.pdf")
    else:
        pyplot.show()
