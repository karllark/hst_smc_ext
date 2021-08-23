# convert Ed Fitzpatrick's IDL save format to the standard
# measure_extincion fits format
import numpy as np
from scipy.io.idl import readsav

import astropy.units as u
from astropy.table import Table

from measure_extinction.extdata import ExtData as ExtDataStock


class ExtData(ExtDataStock):
    """
    Expand the stock ExtData to have a function to read in Ed Fitzpatrick's
    idl save files.
    """

    def read_ext_data_idlsave(self, ext_filename, ebv):
        """
        Read the calculated extinction curve from an idl save file
        In the process, convert the elvebv curves to elv (more fundemental measurement)
        """
        self.type = "elx"
        self.type_rel_band = "V"

        spec_dict = readsav(ext_filename)

        (indxs,) = np.where(np.isfinite(spec_dict["uvcurve"]))
        ngood = len(indxs)
        self.waves["IUE"] = spec_dict["uvwaves"][indxs] * u.Angstrom
        self.exts["IUE"] = spec_dict["uvcurve"][indxs] * ebv
        self.npts["IUE"] = np.full((ngood), 1)
        self.uncs["IUE"] = np.full((ngood), 0.0)

        nphot = len(spec_dict["opirwaves"])
        self.waves["BAND"] = spec_dict["opirwaves"] * u.Angstrom
        self.exts["BAND"] = spec_dict["opircurve"] * ebv
        self.npts["BAND"] = np.full((nphot), 1)
        self.uncs["BAND"] = np.full((nphot), 0.0)
        self.names["BAND"] = ["JohnU", "JohnB", "JohnV", "JohnJ", "JohnH", "JohnK"]

        # model is given as x instead of wavelength
        (indxs,) = np.where(spec_dict["xcurv"] > 0.0)
        self.model["waves"] = np.flip(1.0 / spec_dict["xcurv"][indxs]) * u.micron
        self.model["exts"] = np.flip(spec_dict["bestcurv"][indxs]) * ebv
        nmod = len(self.model["waves"])
        self.model["residuals"] = np.zeros((nmod))
        self.model["type"] = "EdFitzpatrick"


if __name__ == "__main__":

    datapath = "/home/kgordon/Hubble/SMCExt/SENDTOKARL_FINAL/"

    # read in the spreadsheet data
    sdata = Table.read("data/edf_fits_param_edit.csv")

    # read and save the extinction curves
    filelist = "data/smc_stars_all_ed_to_extstar.dat"
    f = open(filelist, "r")
    file_lines = list(f)
    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            names = line.split()
            name = names[0]

            # E(B-V)
            mindx = name.upper() == sdata["name"]
            ebv = float(sdata["ebv"][mindx].data[0])
            ebv_unc = float(sdata["ebv_sig"][mindx].data[0])
            rv = float(sdata["rv"][mindx].data[0])
            rv_unc = float(sdata["rv_sig"][mindx].data[0])
            loghi = float(sdata["logN(HI)"][mindx].data[0])
            loghi_unc = float(sdata["logN(HI)_sig"][mindx].data[0])

            # compute the av
            av = rv * ebv
            if (rv_unc / rv) ** 2 > (ebv_unc / ebv) ** 2:
                av_unc = av * np.sqrt((rv_unc / rv) ** 2 + (ebv_unc / ebv) ** 2)
            else:
                av_unc = av * ebv_unc / ebv

            # read
            text = ExtData()
            text.read_ext_data_idlsave(f"{datapath}{name.upper()}_EXTCURVE.save", ebv)

            # save
            text.columns = {
                "EBV": (ebv, ebv_unc),
                "RV": (rv, rv_unc),
                "LOGHI": (loghi, loghi_unc),
            }
            extcols = {"loghi": loghi, "loghi_unc": loghi_unc}
            savefile = f"fits/{names[1]}_ext.fits"
            text.save(savefile)

            # check the file can be read
            testext = ExtData(filename=savefile)
