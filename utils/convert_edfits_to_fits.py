# convert Ed Fitzpatrick's FITS save format to the standard
# measure_extincion fits format
import numpy as np

import astropy.units as u
from astropy.table import Table
from astropy.io import fits

import matplotlib.pyplot as plt

from measure_extinction.stardata import StarData
from measure_extinction.utils.helpers import get_full_starfile
from measure_extinction.extdata import ExtData as ExtDataStock


class ExtData(ExtDataStock):
    """
    Expand the stock ExtData to have a function to read in Ed Fitzpatrick's
    idl save files.
    """

    def read_ext_data_edfits(self, ext_filename):
        """
        Read the calculated extinction curve from a fits save file
        In the process, convert the elvebv curves to elv (more fundamental measurement)
        """
        self.type = "elx"
        self.type_rel_band = "V"

        # parameters
        # hard coded from eds "FIt_*.out" files
        if "MR12-9" in ext_filename:
            ebv = 0.133
            ebv_unc = 0.003
            loghi = 21.918
            loghi_unc = 0.062
        elif "MR12-10" in ext_filename:
            ebv = 0.249
            ebv_unc = 0.003
            loghi = 22.558
            loghi_unc = 0.192
        elif "MR12-11" in ext_filename:
            ebv = 0.309
            ebv_unc = 0.002
            loghi = 22.003
            loghi_unc = 0.109

        # save
        self.columns = {
            "EBV": (ebv, ebv_unc),
            "LOGHI": (loghi, loghi_unc),
        }

        # UV curve
        edata = fits.getdata(f"{ext_filename}_SPECTROPHOT.fits")
        gvals = np.isfinite(edata[5, :])
        ngood = np.sum(gvals)
        self.waves["IUE"] = edata[0, gvals] * u.Angstrom
        self.exts["IUE"] = edata[5, gvals] * ebv
        self.npts["IUE"] = np.full((ngood), 1)
        self.uncs["IUE"] = np.full((ngood), 0.0)

        # phot curve
        edata = fits.getdata(f"{ext_filename}_FILTERPHOT.fits")
        gvals = np.isfinite(edata[5, :])
        ngood = np.sum(gvals)
        self.waves["BAND"] = edata[0, gvals] * u.Angstrom
        self.exts["BAND"] = edata[5, gvals] * ebv
        self.npts["BAND"] = np.full((ngood), 1)
        self.uncs["BAND"] = np.full((ngood), 0.0)
        self.names["BAND"] = ["WFC3_F225W", "WFC3_F275W", "WFC3_F336W",
                              "ACS_F475W", "ACS_F550M", "ACS_F814W",
                              "WFC3_F110W", "WFC3_F160W"]


if __name__ == "__main__":

    datapath = "/home/kgordon/Python/hst_smc_ext/MR12_ext/"

    # read and save the extinction curves
    filelist = "data/smc_stars_ed_fits_to_extstar.dat"
    f = open(filelist, "r")
    file_lines = list(f)
    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            names = line.split()
            name1 = names[0]
            name2 = names[1]
            print(name1)

            # read
            text = ExtData()
            text.read_ext_data_edfits(f"{datapath}{name1.upper()}")

            # compute the av
            # add estimate uncertainties (not provided by Ed)
            # get the star file to obtain the observed JHK uncs
            fstarname, file_path = get_full_starfile(name2)
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
            text.calc_RV()
            print(text.columns)

            savefile = f"fits/{name2}_ext.fits"
            text.save(savefile)

            # check the file can be read
            testext = ExtData(filename=savefile)
