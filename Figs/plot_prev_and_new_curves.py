import argparse
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

from measure_extinction.extdata import ExtData
from dust_extinction.parameter_averages import G23
from dust_extinction.shapes import FM90

from plot_extinction import mask_bad


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    fontsize = 14

    font = {"size": fontsize}

    plt.rc("font", **font)

    plt.rc("lines", linewidth=1)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("xtick.minor", width=2)
    plt.rc("ytick.major", width=2)
    plt.rc("ytick.minor", width=2)

    figsize = (8.0, 8.0)
    fig, ax = plt.subplots(figsize=figsize)

    g03 = ["prev/azv456_azv70_ext.fits",
           "prev/azv214_azv380_ext.fits",
           "prev/azv18_azv462_ext.fits",
           "prev/azv23_azv404_ext.fits",
           "prev/azv398_azv289_ext.fits"]
    
    g24 = ["fits/azv456_ext_forecor.fits",
           "fits/azv214_ext_forecor.fits",
           "fits/azv18_ext_forecor.fits",
           "fits/azv23_ext_forecor.fits",
           "fits/azv398_ext_forecor.fits"]
    
    g03_labels = ["AzV 456", "AzV 214", "AzV 18", "AzV 23", "AzV 398"]

    offval = 12.0
    for k, cg03 in enumerate(g03):
        pext = ExtData(filename=cg03)
        pext = mask_bad(pext)
        if "IUE" in pext.waves.keys():
            legend_key = "IUE"
        else:
            legend_key = "STIS"

        if k % 2:
            ccol = "b"
        else:
            ccol = "g"

        pext2 = ExtData(filename=g24[k])
        pext2 = mask_bad(pext2)

        pext2.trans_elv_elvebv()
        pext2.plot(ax, 
                   yoffset=k * offval, alpha=0.8, wavenum=True, color=ccol)

        pext.trans_elv_elvebv()
        pext.plot(ax,
                  yoffset=k * offval, 
                  alpha=0.5, 
                  wavenum=True, 
                  color="k",
                  annotate_text=g03_labels[k],
                  annotate_key="BAND",
                  annotate_wave_range=np.array([0.3, 1.0]) * u.micron,
                  annotate_rotation=20.0,
                  annotate_yoffset=1.5)

    # now the MR12 curves
    x = np.arange(3.5,11.0,0.1)/u.micron

    mr12 = [# "MR12-08", 
            "MR12-09", "MR12-10", "MR12-11"]
    fm_params = [# [-2.42, 1.58, 0.44, 4.14, 4.69, 0.53],
                 [0.13, 0.94, 0.32, 5.80, 4.78, 0.68],
                 [-1.51, 1.33, 0.04, 0.23, 4.65, 0.66],
                 [-1.80, 1.44, 1.76, 9.26, 4.76, 0.85]]
    # g24mr12 = [# None,
    #            "fits/mr12-star09_ext_forecor.fits",
    #            "fits/mr12-star10_ext_forecor.fits",
    #            "fits/mr12-star11_ext_forecor.fits"]
    g24mr12 = [# None,
               "fits/mr12-star09_ext.fits",
               "fits/mr12-star10_ext.fits",
               "fits/mr12-star11_ext.fits"]

    for m, cmr12 in enumerate(mr12):
        ext_model = FM90()
        ext_model.parameters = fm_params[m]
        gvals = x >= 6.0 / u.micron
        ax.plot(x[gvals], ext_model(x[gvals]) + (k+m+1) * offval, color="k", alpha=0.5, lw=2, ls="dashed")
        gvals = x < 6.0 / u.micron
        ax.plot(x[gvals], ext_model(x[gvals]) + (k+m+1) * offval, color="k", alpha=0.5, lw=2)        

        if m % 2:
            ccol = "g"
        else:
            ccol = "b"

        if g24mr12[m] is not None:
            mext = ExtData(filename=g24mr12[m])
            mext.rebin_constres("IUE", [0.1, 0.34] * u.micron, 250.)
            mext = mask_bad(mext)
            mext.trans_elv_elvebv()
            mext.plot(ax,
                      yoffset=(k+m+1) * offval, 
                      alpha=0.8, 
                      wavenum=True, 
                      color=ccol,
                      annotate_text=cmr12,
                      annotate_key="BAND",
                      annotate_wave_range=np.array([0.3, 1.0]) * u.micron,
                      annotate_rotation=20.0,
                      annotate_yoffset=1.5)

    ax.set_xlim(0.0, 9.0)
    ax.set_ylim(-5, 120)
    ax.set_ylabel(r"E($\lambda$ - V)/E(B - V) + const")

    # for 2nd x-axis with lambda values
    axis_xs = np.array([0.12, 0.15, 0.2, 0.3, 0.5, 1.0, 3.0])
    new_ticks = 1 / axis_xs
    new_ticks_labels = ["%.2f" % z for z in axis_xs]
    tax = ax.twiny()
    tax.set_xlim(ax.get_xlim())
    tax.set_xticks(new_ticks)
    tax.set_xticklabels(new_ticks_labels, fontsize=0.8*fontsize)
    tax.set_xlabel(r"$\lambda$ [$\mu$m]")

    # ax.legend()

    fig.tight_layout()

    save_str = f"Figs/prev_smc_curves"
    if args.png:
        fig.savefig(f"{save_str}.png")
    elif args.pdf:
        fig.savefig(f"{save_str}.pdf")
    else:
        plt.show()