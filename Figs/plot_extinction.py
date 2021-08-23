import numpy as np
import matplotlib.pyplot as pyplot
import argparse
import matplotlib

from scipy.io.idl import readsav
import astropy.units as u

from measure_extinction.extdata import ExtData as ExtDataStock


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


def plot_ext_stack(filelist, ax, idlsave=False, locpath="./", fontsize=14):

    f = open(filelist, "r")
    file_lines = list(f)
    starnames = []
    extdatas = []
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
                text.read(locpath + starnames[-1] + "_ext.fits")
            extdatas.append(text)

    ann_wave_range = [15.0, 18.0]
    col_vals = ["b", "g", "c"]
    # lin_vals = ["--", ":", "-."]
    n_cols = len(col_vals)

    ann_wave_range = 1.0 / np.array([0.3, 0.25]) / u.micron
    ann_wave_range = np.array([0.3, 3.]) / u.micron
    n_stars = len(starnames)
    offset_val = 5.0

    for k, cdata in enumerate(extdatas):

        cdata.trans_elv_elvebv()

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
            model=True,
        )

    ax.set_yscale("linear")

    ax.set_xlabel(r"$1/\lambda$ [$\mu m^{-1}$]")

    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax.set_xlim(1.0 / 2.5, 10.5)
    ax.set_ylim(-5.0, 20.0 + offset_val * n_stars)


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
