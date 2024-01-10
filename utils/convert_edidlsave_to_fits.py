# convert Ed Fitzpatrick's IDL save format to the standard
# measure_extincion fits format
import numpy as np
from scipy.io import readsav

import astropy.units as u
from astropy.table import Table

from measure_extinction.stardata import StarData
from measure_extinction.utils.helpers import get_full_starfile
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
            name1 = names[0]  # name of IDL save file
            name2 = names[1]  # name in CSV file
            name3 = names[2]  # homogenized name for paper

            # E(B-V)
            mindx = name2 == sdata["name"]
            ebv = float(sdata["ebv"][mindx].data[0])
            ebv_unc = float(sdata["ebv_sig"][mindx].data[0])
            rv = float(sdata["rv"][mindx].data[0])
            rv_unc = float(sdata["rv_sig"][mindx].data[0])
            loghi = float(sdata["logN(HI)"][mindx].data[0])
            loghi_unc = float(sdata["logN(HI)_sig"][mindx].data[0])

            # read
            text = ExtData()
            text.read_ext_data_idlsave(f"{datapath}{name1.upper()}_EXTCURVE.save", ebv)

            # compute the av
            # add estimate uncertainties (not provided by Ed)
            # get the star file to obtain the observed JHK uncs
            fstarname, file_path = get_full_starfile(name3)
            starobs = StarData(fstarname, path=file_path)
            
            cmags = ["J", "H", "K"]
            cwaves = np.array([1.25, 1.6, 2.2]) * u.micron
            for cmag, cwave in zip(cmags, cwaves):
                exv_unc = np.sqrt(starobs.data["BAND"].bands[cmag][1]**2
                                + starobs.data["BAND"].bands["V"][1]**2)
                dwaves = np.absolute(text.waves["BAND"] - cwave)
                kindx = dwaves.argmin()
                if dwaves[kindx] < 0.1 * u.micron:
                    text.uncs["BAND"][kindx] = exv_unc

            text.calc_AV_JHK()
            text.columns["EBV"] = (ebv, ebv_unc)
            text.calc_RV()

            print("RV Ed, Karl", name3, rv, rv_unc, text.columns["RV"])

            # save
            text.columns["LOGHI"] = (loghi, loghi_unc)
            savefile = f"fits/{name3}_ext.fits"
            text.save(savefile)

            # check the file can be read
            testext = ExtData(filename=savefile)
