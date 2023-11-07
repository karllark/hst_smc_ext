import numpy as np
import matplotlib.pyplot as pyplot
import argparse
import matplotlib

from scipy.io.idl import readsav
import astropy.units as u

from measure_extinction.extdata import ExtData as ExtDataStock


def mask_bad(cdata):
    # mask bad data (lines)
    mask = [(8.4, 8.05),
            (7.55, 7.45),
            (7.75, 7.65),
            (7.95, 7.9),
            (6.5, 6.4),
            (4.22, 4.17),
            (4.28, 4.255),
            (3.51, 3.49),
            (3.59, 3.56),
            (3.87, 3.82)]
    mask = 1.0 / np.array(mask)
    for region in mask:
        for src in cdata.waves.keys():
            cdata.npts[src][
                (cdata.waves[src].value >= region[0])
                & (cdata.waves[src].value <= region[1])
            ] = 0
    return cdata

class ExtData(ExtDataStock):
    """
    Expand the stock ExtData to have a function to read in Ed Fitzpatrick's
    idl save files.
    """

    def read_ext_data_idlsave(self, ext_filename):
        """
        Read the calculated extinction curve from an idl save file
        """
        spec_dict = readsav(ext_filename)

        (indxs,) = np.where(np.isfinite(spec_dict["uvcurve"]))
        ngood = len(indxs)
        self.waves["STIS"] = spec_dict["uvwaves"][indxs] * u.Angstrom
        self.exts["STIS"] = spec_dict["uvcurve"][indxs]
        self.npts["STIS"] = np.full((ngood), 1)
        self.uncs["STIS"] = np.full((ngood), 0.0)

        nphot = len(spec_dict["opirwaves"])
        self.waves["BAND"] = spec_dict["opirwaves"] * u.Angstrom
        self.exts["BAND"] = spec_dict["opircurve"]
        self.npts["BAND"] = np.full((nphot), 1)
        self.uncs["BAND"] = np.full((nphot), 0.0)

        # model is given as x instead of wavelength
        (indxs,) = np.where(spec_dict["xcurv"] > 0.0)
        self.waves["MODEL"] = np.flip(1.0 / spec_dict["xcurv"][indxs]) * u.micron
        self.exts["MODEL"] = np.flip(spec_dict["bestcurv"][indxs])
        nmod = len(self.waves["MODEL"])
        self.npts["MODEL"] = np.full((nmod), 1)
        self.uncs["MODEL"] = np.full((nmod), 0.0)


def plot_ext_stack(
    filelist, ax, idlsave=False, locpath="./", fontsize=14, forecor=False, adjusted=False,
):

    f = open(filelist, "r")
    file_lines = list(f)
    starnames = []
    extdatas = []
    spslopes = []
    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            starnames.append(name)
            text = ExtData()
            if idlsave:
                text.read_ext_data_idlsave(
                    locpath
                    + starnames[-1].upper()
                    + "_EXTCURVE.save"
                    #                    locpath + "FITSPEC_" + starnames[-1].upper() + "_ATLAS.save"
                )
            else:
                filebase = locpath + starnames[-1] + "_ext"
                if forecor:
                    filebase = f"{filebase}_forecor"
                if adjusted:
                    filebase = f"{filebase}_adjusted"
                print(filebase)
                text.read(f"{filebase}.fits")
            extdatas.append(text)

            if "IUE" not in text.waves.keys():
                bkey = "STIS"
            else:
                bkey = "IUE"

            gvals1 = np.logical_and(
                text.waves[bkey] > 0.14 * u.micron,
                text.waves[bkey] < 0.20 * u.micron,
            )
            gvals1 = np.logical_and(gvals1, text.npts[bkey] > 0)
            gvals2 = np.logical_and(
                text.waves[bkey] > 0.22 * u.micron,
                text.waves[bkey] < 0.30 * u.micron,
            )
            gvals2 = np.logical_and(gvals2, text.npts[bkey] > 0)
            spslopes.append(
                np.median(text.exts[bkey][gvals1]) - np.median(text.exts[bkey][gvals2])
            )

    slpsort = np.argsort(spslopes)

    ann_wave_range = [15.0, 18.0]
    col_vals = ["b", "g", "c"]
    # lin_vals = ["--", ":", "-."]
    n_cols = len(col_vals)

    ann_wave_range = 1.0 / np.array([0.3, 0.25]) / u.micron
    ann_wave_range = np.array([0.3, 3.0]) / u.micron
    n_stars = len(starnames)
    offset_val = 5.0

    # for k, cdata in enumerate(extdatas):
    for j in range(len(extdatas)):
        # k = slpsort[j]
        # print(k)
        k = j
        cdata = extdatas[k]

        cdata.trans_elv_elvebv()

        if "mr12" in starnames[k]:
            rebinfac = 5
        else:
            rebinfac = None

        # mask bad data (lines)
        mask = [(8.4, 8.05),
                (7.55, 7.45),
                (7.75, 7.65),
                (7.95, 7.9),
                (6.5, 6.4),
                (4.22, 4.17),
                (4.28, 4.255),
                (3.51, 3.49),
                (3.59, 3.56),
                (3.87, 3.82)]
        mask = 1.0 / np.array(mask)
        for region in mask:
            for src in cdata.waves.keys():
                cdata.npts[src][
                    (cdata.waves[src].value >= region[0])
                    & (cdata.waves[src].value <= region[1])
                ] = 0

        # plot the extinction curves
        cdata.plot(
            ax,
            color=col_vals[k % n_cols],
            yoffset=k * offset_val,
            alpha=0.5,
            fontsize=0.75 * fontsize,
            wavenum=True,
            annotate_text=starnames[k],
            annotate_key="BAND",
            annotate_wave_range=ann_wave_range,
            annotate_rotation=10.0,
            annotate_yoffset=1.5,
            model=False,
            rebin_fac=rebinfac,
        )

    ax.set_yscale("linear")

    # ax.set_xlim(1.0 / 2.5, 10.5)
    ax.set_ylim(-5.0, 30.0 + offset_val * n_stars)

    ax.set_xlabel(r"$1/\lambda$ [$\mu m^{-1}$]", fontsize=fontsize)

    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")

    # for 2nd x-axis with lambda values
    axis_xs = np.array([0.12, 0.15, 0.2, 0.3, 0.5, 1.0])
    new_ticks = 1 / axis_xs
    new_ticks_labels = ["%.2f" % z for z in axis_xs]
    tax = ax.twiny()
    tax.set_xlim(ax.get_xlim())
    tax.set_xticks(new_ticks)
    tax.set_xticklabels(new_ticks_labels, fontsize=0.8*fontsize)
    tax.set_xlabel(r"$\lambda$ [$\mu$m]")

    # ax.spines["right"].set_visible(False)
    # ax.spines["top"].set_visible(False)


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of stars to plot")
    parser.add_argument(
        "-i",
        "--idlsave",
        help="extinction data in an IDL save file",
        action="store_true",
    )
    parser.add_argument(
        "-p", "--png", help="save figure as a png file", action="store_true"
    )
    parser.add_argument(
        "-e", "--eps", help="save figure as an eps file", action="store_true"
    )
    parser.add_argument("--path", help="path to filelist and extinction curves")
    args = parser.parse_args()

    if args.path:
        locpath = args.path + "/"
    else:
        locpath = ""

    fontsize = 18

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    fig, ax = pyplot.subplots(nrows=1, ncols=1, figsize=(10, 12))

    plot_ext_stack(args.filelist, ax, idlsave=args.idlsave, locpath=locpath)

    ax.set_xlabel(r"$1/\lambda$ [$\mu m^{-1}$]", fontsize=1.3 * fontsize)
    ax.set_ylabel(r"$E(\lambda - V)/E(B - V)$ + offset", fontsize=1.3 * fontsize)

    pyplot.show()
