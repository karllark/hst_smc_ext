# convert Ed Fitzpatrick's IDL save format to the standard
# measure_extincion fits format
import numpy as np

import astropy.units as u
from astropy.table import Table

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
        rv = ptab["R(55)"][0]
        rv_unc = 0.0
        loghi = ptab["log(NHI)"][0]
        loghi_unc = 0.0

        # compute the av
        av = rv * ebv
        if (rv_unc / rv) ** 2 > (ebv_unc / ebv) ** 2:
            av_unc = av * np.sqrt((rv_unc / rv) ** 2 + (ebv_unc / ebv) ** 2)
        else:
            av_unc = av * ebv_unc / ebv

        # save
        self.columns = {
            "EBV": (ebv, ebv_unc),
            "RV": (rv, rv_unc),
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
            name = names[0]
            print(name)

            # read
            text = ExtData()
            text.read_ext_data_edtxt(f"{datapath}{name.upper()}")

            savefile = f"fits/{names[1]}_ext.fits"
            text.save(savefile)

            # check the file can be read
            testext = ExtData(filename=savefile)
