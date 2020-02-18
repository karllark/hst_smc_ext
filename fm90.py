# routine to compute the FM90 UV extinction curve (6 parameters)

import numpy as np

def fm90(params, x):

    x2 = x**2
    elvebv = (params[0] + params[1]*x +  # linear term
              params[2]*(x2/((x2 - params[4]**2)**2 + x2*(params[5]**2))))  # 2175 A bump term

    fnuv_indxs = np.where(x >= 5.9)
    if len(fnuv_indxs) > 0:
        y = x[fnuv_indxs] - 5.9
        elvebv[fnuv_indxs] += params[3]*(0.5392*(y**2) + 0.05644*(y**3))  # FUV curvature term

    return elvebv
