# convert Ed Fitzpatrick's IDL save format to the standard
# measure_extincion fits format
import numpy as np

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

    def read_ext_data_edtxt(self, ext_filename):
        """
        Read the calculated extinction curve from an idl save file
        In the process, convert the elvebv curves to elv (more fundemental measurement)
        """
        self.type = "elx"
        self.type_rel_band = "V"

        # parameters
        ptab = Table.read(f"{ext_filename}_parms.txt", format="ascii.commented_header", header_start=1)
        ebv = ptab["E(44-55)"][0]
        ebv_unc = 0.0
        loghi = ptab["log(NHI)"][0]
        loghi_unc = 0.0

        # save
        self.columns = {
            "EBV": (ebv, ebv_unc),
            "LOGHI": (loghi, loghi_unc),
        }
        # extcols = {"loghi": loghi, "loghi_unc": loghi_unc}

        # UV curve
        spec_dict = Table.read(f"{ext_filename}_uvcurv.txt", format="ascii.commented_header", header_start=1)

        (indxs,) = np.where(np.isfinite(spec_dict["K(LAMBDA-55)"]))
        ngood = len(indxs)
        self.waves["IUE"] = spec_dict["LAMBDA"][indxs] * u.Angstrom
        self.exts["IUE"] = spec_dict["K(LAMBDA-55)"][indxs] * ebv
        self.npts["IUE"] = np.full((ngood), 1)
        self.uncs["IUE"] = np.full((ngood), 0.0)

        # UV curve
        spec_dict = Table.read(f"{ext_filename}_opircurv.txt", format="ascii.commented_header", header_start=1)

        nphot = len(spec_dict["LAMBDA(EFF)"])
        self.waves["BAND"] = spec_dict["LAMBDA(EFF)"] * u.Angstrom
        self.exts["BAND"] = spec_dict["K(LAMBDA-55)"] * ebv
        self.npts["BAND"] = np.full((nphot), 1)
        self.uncs["BAND"] = np.full((nphot), 0.0)
        self.names["BAND"] = ["WFC3_F336W", "ACS_F475W", "ACS_F550M", "ACS_F814W", "WFC3_F110W", "WFC3_F160W"]


if __name__ == "__main__":

    datapath = "/home/kgordon/Python/hst_smc_ext/MR12_ext/"

    # read and save the extinction curves
    filelist = "data/smc_stars_ed_txt_to_extstar.dat"
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
            text.read_ext_data_edtxt(f"{datapath}{name1.upper()}")

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
