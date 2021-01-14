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
to read the data and calculate extinction curves.

Figures
-------

1. ??

2. UV spectra of all program stars: Figs/plot_uv_mspec.py data/smc_stars_all.dat

3. from Ralph Bohlin

4. Figs/plot_smc_3panel.py

5. Figs/model_contamination.py

Tables
------

1.
