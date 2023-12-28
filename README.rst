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

4. UV extinction curves for set A: Figs/plot_smc_2panel.py --both

5. Comparison of new to previous G03, MR12 UV curves:
   Figs/plot_prev_and_new_curves.py

6. Group type versus position in SMC: Figs/plot_positions_mips24.py
   Plot the bump/nobump/flat/lowebv samples over the MIPS 24um image 

7. FM90 parameters versus E(B-V) etc: TBR
   all these done with extinction_ensemble_parameters repository

8. Bump versus q_PAH: Figs/plot_bump_vs_qpah.py

9. SMC Average (& Bump Average) curves versus MW/LMC: Figs/plot_smcave.py
   add `--bumps` to get the average of the sightlines with bumps

Tables
------

1. Manually

2. Manually - from Ralph

3. MW foreground HI and E(B-V) & q_PAH: From Figs/plot_samp_prop.py
   output as data/ancillary_samp_properties.tex
   qpah numbers table created when running Figs/plot_positions_mips24.py

4. Stellar fit parameters: utils/create_stellar_param_table.py

5. Sightline properties: utils/create_param_table.py
   E(B-V), R(V), N(HI), FM90 parameters
   
6. SMC Average/Bump Average versus wavelength: TBC  

Appendix Tables
---------------

1. Manually - from Ed

2. Manually - from Ed

