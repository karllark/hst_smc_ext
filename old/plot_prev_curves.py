#!/usr/bin/env python
#
# Program to plot a list of extinction curves
#
# written Apr 2015 by Karl Gordon (kgordon@stsci.edu)
# based strongly on IDL programs created over the previous 10 years
#
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib
import argparse

from extdata import *
from fm90 import fm90
from f99 import f99

if __name__ == "__main__":

    font = {'family' : 'normal',
            'size'   : 15}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of stars to plot")
    parser.add_argument("-d", "--data", help="plot data instead of FM90 fit",
                        action="store_true")
    parser.add_argument("-p", "--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("-e", "--eps", help="save figure as an eps file",
                        action="store_true")
    parser.add_argument("--path", help="path to filelist and extinction curves")
    args = parser.parse_args()

    if args.path:
        locpath = args.path + '/'
    else:
        locpath = ''

    f = open(args.filelist, 'r')
    file_lines = list(f)
    starnames = []
    extdata = []
    for line in file_lines:
        if (string.find(line,'#') != 0) & (len(line) > 0):
            name = line.rstrip()
            starnames.append(name)
            text = ExtData()
            text.read_ext_data(starnames[-1]+'_ext.fits')
            extdata.append(text)
            
    fig = pyplot.figure(figsize=(12.,12.))
    ax = fig.add_subplot(1,1,1)

    off1 = 10.
    off2 = 20.
    off3 = 30.

    mw31_x = np.arange(0.1,9.0,0.02)
    mw31_ext = f99(3.1, mw31_x)
    mw31_ext = (mw31_ext - 1.0)*3.1
    ax.plot(mw31_x, mw31_ext, 'k--')
    ax.plot(mw31_x, mw31_ext+off1, 'k--')
    ax.plot(mw31_x, mw31_ext+off2, 'k--')
    ax.plot(mw31_x, mw31_ext+off3, 'k--')

    norm_wave_range = [0.3, 1.0]
    ann_wave_range = 1.0/np.array([8.65,8.45])
    col_vals = ['b','g','m','c','y']
    offset_val = 5.0

    offset_vals = np.concatenate((np.full(4, 0.0),[off1],np.full(11, off2),np.full(4, off3)))

    for i in range(len(starnames)):
        k = i
        offset_val = offset_vals[i]

        vindxs = np.argsort(abs(extdata[k].ext_waves['BANDS'] - 0.545))
        vk = vindxs[0]
        bindxs = np.argsort(abs(extdata[k].ext_waves['BANDS'] - 0.438))
        bk = bindxs[0]

        ebvebv = extdata[k].ext_curve['BANDS'][bk] - extdata[k].ext_curve['BANDS'][vk]
        print(starnames[i], extdata[k].ext_curve['BANDS'][bk] - extdata[k].ext_curve['BANDS'][vk], extdata[k].ebv)
        norm_val = 1.0/ebvebv

        if k < 4:
            ck = 0
        elif k < 5:
            ck = 1
        else:
            ck = 2

        if args.data:
            ax.plot(1.0/extdata[k].ext_waves['BANDS'],extdata[k].ext_curve['BANDS']*norm_val+offset_val,col_vals[ck]+'o')
            ax.plot(1.0/extdata[k].ext_waves['STIS'],extdata[k].ext_curve['STIS']*norm_val+offset_val,col_vals[ck]+'-')
        else:
            ax.plot(1.0/extdata[k].ext_waves['BANDS'],extdata[k].ext_curve['BANDS']*norm_val+offset_val,col_vals[ck]+'o')

            c2 = extdata[k].fm90['C2'][0]
            c1 = 2.09 - 2.84*c2  # from FM07
            fm_params = [c1,c2,extdata[k].fm90['C3'][0],extdata[k].fm90['C4'][0],
                         extdata[k].fm90['x0'][0],extdata[k].fm90['gamma'][0]]
    
            elvebv = fm90(fm_params, 1.0/extdata[k].ext_waves['STIS'])
            ax.plot(1.0/extdata[k].ext_waves['STIS'],elvebv*norm_val+offset_val,col_vals[ck]+'-')

    # plot M...
    offset_val = 30

    gindxs, = np.where(1.0/extdata[0].ext_waves['STIS'] < 6.0)
    bindxs, = np.where(1.0/extdata[0].ext_waves['STIS'] >= 6.0)

    fm_params = [-2.42, 1.58, 0.44, 4.14, 4.69, 0.53]
    elvebv = fm90(fm_params, 1.0/extdata[0].ext_waves['STIS'])
    ax.plot(1.0/extdata[0].ext_waves['STIS'][gindxs],elvebv[gindxs]+offset_val,col_vals[3]+'-')
    ax.plot(1.0/extdata[0].ext_waves['STIS'][bindxs],elvebv[bindxs]+offset_val,col_vals[3]+'--')

    fm_params = [0.13, 0.94, 0.32, 5.80, 4.78, 0.68]
    elvebv = fm90(fm_params, 1.0/extdata[0].ext_waves['STIS'])
    ax.plot(1.0/extdata[0].ext_waves['STIS'][gindxs],elvebv[gindxs]+offset_val,col_vals[3]+'-')
    ax.plot(1.0/extdata[0].ext_waves['STIS'][bindxs],elvebv[bindxs]+offset_val,col_vals[3]+'--')

    fm_params = [-1.51, 1.33, 0.04, 0.23, 4.65, 0.66]
    elvebv = fm90(fm_params, 1.0/extdata[0].ext_waves['STIS'])
    ax.plot(1.0/extdata[0].ext_waves['STIS'][gindxs],elvebv[gindxs]+offset_val,col_vals[3]+'-')
    ax.plot(1.0/extdata[0].ext_waves['STIS'][bindxs],elvebv[bindxs]+offset_val,col_vals[3]+'--')

    fm_params = [-1.80, 1.44, 1.76, 9.26, 4.76, 0.85]
    elvebv = fm90(fm_params, 1.0/extdata[0].ext_waves['STIS'])
    ax.plot(1.0/extdata[0].ext_waves['STIS'][gindxs],elvebv[gindxs]+offset_val,col_vals[3]+'-')
    ax.plot(1.0/extdata[0].ext_waves['STIS'][bindxs],elvebv[bindxs]+offset_val,col_vals[3]+'--')


    #ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_ylim(-5.,60.)
    ax.set_xlim(1.0/3.0,1.0/0.1)

    fontsize = 22
    ax.set_xlabel('$1/\lambda [\mu m^{-1}]$', fontsize=fontsize)
    ax.set_ylabel('$E(\lambda - V)/E(B - V) + \mathrm{offset}$', fontsize=fontsize)

    ax.text(6., 3.0, "$\mathrm{MW}$ $R(V) = 3.1$", color='k')
    #ax.text(7.0, 5., "AzV 456", color=col_vals[0])
    
    tx = 9.5
    dtx = 0.25
    ax.text(tx, 10., "Gordon & Clayton", color=col_vals[0], rotation=270., verticalalignment='center')
    ax.text(tx-dtx, 10., "(1998)", color=col_vals[0], rotation=270., verticalalignment='center')

    ax.text(tx, 22., "Gordon et al.", color=col_vals[1], rotation=270., verticalalignment='center')
    ax.text(tx-dtx, 22., "(2003)", color=col_vals[1], rotation=270., verticalalignment='center')

    ax.text(tx, 33., "Gordon et al.", color=col_vals[2], rotation=270., verticalalignment='center')
    ax.text(tx-dtx, 33., "(2015)", color=col_vals[2], rotation=270., verticalalignment='center')

    ax.text(tx, 50., "Maiz Apellaniz & Rubio", color=col_vals[3], rotation=270., verticalalignment='center')
    ax.text(tx-dtx, 50., "(2012)", color=col_vals[3], rotation=270., verticalalignment='center')

    pyplot.tight_layout()

    save_str = '_ext_hst23'
    if args.png:
        fig.savefig(string.replace(args.filelist,'.dat',save_str+'.png'))
    elif args.eps:
        fig.savefig(string.replace(args.filelist,'.dat',save_str+'.eps'))
    else:
        pyplot.show()

