# routine to compute the F99 R(V) dependent extinction curve
#  with the option to change the UV FM90 parameters

import numpy as np
from scipy import interpolate
import fm90

def f99(Rv, x, x0=4.596, gamma=0.99, c3=3.23, c4=0.41, c2=None):

    # default FM90 parameters set in function definition and here
    if c2 == None:
        c2 = -0.824 + 4.717/Rv

    c1 = 2.09 - 2.84*c2  # updated for FM07

    # initialize
    curve = np.zeros(len(x))

    # Compute UV portion of A(lambda)/E(B-V) curve using FM fitting function and 
    # R-dependent coefficients
    x_cutval_uv = 10000.0/2700.0
    x_splineval_uv = 10000.0/np.array([2700.0,2600.0])
    indxs_uv = np.where(x >= x_cutval_uv)
    if len(indxs_uv) > 0:
        xuv = np.concatenate([x_splineval_uv, x[indxs_uv]])
    else:
        xuv = x_splineval_uv

    yuv = fm90.fm90([c1,c2,c3,c4,x0,gamma],xuv)
    yuv = yuv + Rv  # E(l-V)/E(B-V) to A(V)/E(B-V)

    y_splineval_uv = yuv[0:2]              # save spline points

    if len(indxs_uv) > 0:
        curve[indxs_uv] = yuv[2:]   # ingore the spline points

    # Compute optical portion of A(lambda)/E(B-V) curve
    # using cubic spline anchored in UV, optical, and IR
    indxs_opir = np.where(x < x_cutval_uv)

    x_splineval_optir = 10000./np.array([26500.0,12200.0,6000.0,5470.0,4670.0,4110.0])
    x_splineval_optir = np.insert(x_splineval_optir, 0, 0.0)

    y_splineval_ir    = np.array([0.0,0.265,0.829])*Rv/3.1 
    y_splineval_opt   = np.array([-0.426 + 1.0044*Rv,
                                  -0.050 + 1.0016*Rv,
                                  0.701 + 1.016*Rv,
                                  1.208 + 1.0032*Rv - 0.00033*(Rv**2)])

    y_splineval_optir = np.concatenate([y_splineval_ir,y_splineval_opt])

    if len(indxs_opir) > 0:
        spline_x = np.concatenate([x_splineval_optir, x_splineval_uv])
        spline_y = np.concatenate([y_splineval_optir, y_splineval_uv])
        spline_rep = interpolate.splrep(spline_x, spline_y)
        curve[indxs_opir] = interpolate.splev(x[indxs_opir], spline_rep, der=0)

    return curve/Rv
