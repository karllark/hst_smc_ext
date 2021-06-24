Code for a SMC UV Extinction Paper
==================================

Routines for SMC UV dust extinction curve paper.
This paper is Gordon, Fitzpatrick, et al. (in prep).

In Development!
---------------

Active development.
Data still in changing.
Use at your own risk.

Contributors
------------
Karl Gordon

License
-------

This code is licensed under a 3-clause BSD style license (see the
``LICENSE`` file).

Extinction Curves
-----------------

External packages used include `dust_extintion`, `measure_extinction`, and
`extstar_data`.  The basic information on each star is given in `extstar_data`
including (TBD) the STIS spectra.  `dust_extinction` provides the
fitting curves shapes (FM90, AxAvToExv).  `measure_extinction` gives the routines
to read the data and extinction curves.

Extinction curves measured by Ed Fitzpatrick using his custom IDL routines.
These curves were converted to the standard `measure_extinction` format using
utils/convert_edidlsave_to_fits.py.

Generate full list of extinction curves sightlines (in data dir).
`cat smc_stars_reddened_good_highebv.dat smc_stars_reddened_good_lowebv.dat smc_stars_reddened_suspect.dat | sort > smc_stars_all.dat`

Converting the extinction curves measured by Ed F.
`utils/convert_edidlsave_to_fits.py`

Figures
-------

1. ?? (from proposal - revision needed)

2. UV spectra of all program stars: Figs/plot_uv_mspec.py data/smc_stars_all.dat

3. from Ralph Bohlin

4. Figs/plot_smc_3panel.py

5. Figs/model_contamination.py

6. Figure showing the obspair, modelpair, and difference between them.
   This is measureing the MW foreground extinction.
   Figs/plot_old_new.py

7. Figure showing the normalized MW foreground extinctions.
   Fig/plot_mwforeground.py

Tables
------

1.
