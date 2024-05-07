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
including the STIS spectra.  `dust_extinction` provides the
fitting curves shapes (FM90_B3).  `measure_extinction` gives the routines
to read the data and extinction curves.

Extinction curves measured by Ed Fitzpatrick using his custom IDL routines.
These curves were converted to the standard `measure_extinction` format using
utils/convert_edidlsave_to_fits.py.

Foreground correction done with Figs/plot_forecor_ext.py.  Script to do all the
sightlines (including the adjusted versions) is forecor_all.

Figures
-------

1. UV spectra of all program stars: Figs/plot_uv_mspec_byprog.py

2. from Ralph Bohlin

3. Examples of foreground subtraction:
   Figs/plot_forecor_ext.py azv23 0.321  0.064
   Figs/plot_forecor_ext.py azv456 0.267  0.031

4. UV extinction curves: Figs/plot_smc_2panel.py --both

5. Comparison of new to previous G03, MR12 UV curves:
   Figs/plot_prev_and_new_curves.py

6. Group type versus position in SMC: Figs/plot_positions_mips24.py
   Plot the bump/nobump/flat/lowebv samples over the MIPS 24um image 

7. SMC Average (& Bump Average) curves versus MW/LMC: Figs/plot_smcave.py
   add `--bumps` to get the average of the sightlines with bumps

8. Bump versus q_PAH: Figs/plot_bump_vs_qpah.py --inclmc

9. R(V) versus A(V) and N(HI) versus A(V): extinction_ensemble_parameters repository
   plot_param_vs_param.py --datasets gor09 gor03_lmc gor24_smc_nobump gor24_smc_bump gor24_smc_flat --spropsav

10. FM90 parameters versus A(V): extinction_ensemble_parameters repository
   plot_param_vs_param.py --datasets gor09 gor03_lmc gor24_smc_nobump gor24_smc_bump gor24_smc_flat --av

11. FM90 parameters versus each other: extinction_ensemble_parameters repository
   plot_param_vs_param.py --datasets gor09 gor03_lmc gor24_smc_nobump gor24_smc_bump --fm90main --fit

12. A(V) and FM90 versus N(HI)/A(V): extinction_ensemble_parameters repository
     plot_param_vs_param.py --datasets gor09 gor03_lmc gor24_smc_nobump gor24_smc_bump --gdprops

Tables
------

1. Manually

2. Manually - from Ralph

3. MW foreground HI and E(B-V) & q_PAH: From Figs/plot_samp_prop.py
   output as data/ancillary_samp_properties.tex
   qpah numbers table created when running Figs/plot_positions_mips24.py

4. Stellar fit parameters: utils/create_stellar_param_table.py

5. & 6. Sightline properties: utils/create_param_table.py
   E(B-V), R(V), N(HI), FM90 parameters
   
7. SMC Average/Bump Average versus wavelength: utils/create_aveext_table.py
   Uses ASCII files created by Figs/plot_smcave.py

Appendix Tables
---------------

1. UVB photometry: Manually - from Ed

2. griz and JHK_S photometry: Manually - from Ed

3. MR12 HST photometry: utils/create_mr12_photom_table.py

