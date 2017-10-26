"""
desispec.flexure
==================

Functions to measure 'flexure' in the instrument using
the night sky lines

A fair bit of this has been ported/grabbed from PYPIT
"""
from __future__ import print_function, absolute_import, division

import numpy as np
import copy
import pdb
import imp
import yaml
import glob
import math
import time
import os
import sys
import argparse
import locale
from pkg_resources import resource_exists, resource_filename

from astropy.modeling import models, fitting
from astropy.stats import sigma_clip
from astropy.table import Table, Column, vstack
from astropy.io import fits

from desispec.util import set_backend
set_backend()

from matplotlib import pyplot as plt
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages

from desiutil.log import get_logger
from desispec.spectrum import Spectrum
from desispec.bootcalib import find_arc_lines



def flexure_archive():
    """  Load archived sky spectrum
    Returns:
        spec -- Spectrum object
    """
    #   latitude = settings.spect['mosaic']['latitude']
    #   longitude = settings.spect['mosaic']['longitude']
    skyspec_fil = resource_filename('desispec', 'data/sky/paranal_sky.fits')
    hdu = fits.open(skyspec_fil)
    # Instantiate
    spec = Spectrum(hdu[1].data, hdu[0].data, np.ones_like(hdu[0].data))
    #
    return spec


def flex_shift(channel, obj_skyspec, arx_skyspec, dwv_arx=0.05):
    """ Calculate shift between object sky spectrum and archive sky spectrum

    import desimodel.io
    psf = desimodel.io.load_psf('b')
    w  = np.linspace(psf.wmin, psf.wmax, 100)
    wd = psf.wdisp(100, w)
    np.median(wd)


    Parameters
    ----------
    slf
    det
    obj_skyspec
    arx_skyspec
    dwv_arx : float
      Dispersion of Paranal spectrum
      R_Paranal = 31300

    Returns
    -------
    flex_dict
    """
    # DESI resoultion
    wdisp_dict = dict(b=0.628, r=0.578, z=0.731)  # Ang
    # Paranal
    wdisp_paranal =

    # Find the brightest emission lines
    arx_cent, arx_amp = find_arc_lines(arx_skyspec.flux)
    obj_cent, obj_amp = find_arc_lines(obj_skyspec.flux)

    #Keep only 5 brightest amplitude lines (xxx_keep is array of indices within arx_w of the 5 brightest)
    arx_keep = np.argsort(arx_amp)[-5:]
    obj_keep = np.argsort(obj_amp)[-5:]

    #Calculate dispersion (Angstrom per pixel)
    arx_disp = np.append(arx_skyspec.wavelength[1]-arx_skyspec.value[0],
                         arx_skyspec.wavelength[1:]-arx_skyspec.value[:-1])
    obj_disp = np.append(obj_skyspec.wavelength[1]-obj_skyspec.wavelength[0],
                         obj_skyspec.wavelength[1:]-obj_skyspec.wavelength[:-1])

    '''
    #Calculate resolution (lambda/delta lambda_FWHM)..maybe don't need this? can just use sigmas
    arx_idx = (arx_cent+0.5).astype(np.int)[arx_keep]   # The +0.5 is for rounding
    arx_res = arx_skyspec.wavelength[arx_idx]/\
              (arx_disp[arx_idx]*(2*np.sqrt(2*np.log(2)))*arx_wid[arx_w][arx_keep])
    obj_idx = (obj_cent+0.5).astype(np.int)[obj_w][obj_keep]   # The +0.5 is for rounding
    obj_res = obj_skyspec.wavelength.value[obj_idx]/ \
              (obj_disp[obj_idx]*(2*np.sqrt(2*np.log(2)))*obj_wid[obj_w][obj_keep])
    #obj_res = (obj_sky.wavelength.value[0]+(obj_disp*obj_cent[obj_w][obj_keep]))/(
    #    obj_disp*(2*np.sqrt(2*np.log(2)))*obj_wid[obj_w][obj_keep])
    msgs.info("Resolution of Archive={:g} and Observation={:g}".format(
        np.median(arx_res), np.median(obj_res)))
    '''

    #Determine sigma of gaussian for smoothing
    arx_sig2 = (arx_disp[arx_idx]*arx_wid[arx_w][arx_keep])**2.
    obj_sig2 = (obj_disp[obj_idx]*obj_wid[obj_w][obj_keep])**2.

    arx_med_sig2 = np.median(arx_sig2)
    obj_med_sig2 = np.median(obj_sig2)

    if obj_med_sig2 >= arx_med_sig2:
        smooth_sig = np.sqrt(obj_med_sig2-arx_med_sig2)  # Ang
        smooth_sig_pix = smooth_sig / np.median(arx_disp[arx_idx])
        arx_skyspec = arx_skyspec.gauss_smooth(smooth_sig_pix*2*np.sqrt(2*np.log(2)))
    else:
        msgs.warn("Prefer archival sky spectrum to have higher resolution")
        smooth_sig_pix = 0.
        msgs.warn("New Sky has higher resolution than Archive.  Not smoothing")
        #smooth_sig = np.sqrt(arx_med_sig**2-obj_med_sig**2)

    #Determine region of wavelength overlap
    min_wave = max(np.amin(arx_skyspec.wavelength.value), np.amin(obj_skyspec.wavelength.value))
    max_wave = min(np.amax(arx_skyspec.wavelength.value), np.amax(obj_skyspec.wavelength.value))

    #Smooth higher resolution spectrum by smooth_sig (flux is conserved!)
#    if np.median(obj_res) >= np.median(arx_res):
#        msgs.warn("New Sky has higher resolution than Archive.  Not smoothing")
        #obj_sky_newflux = ndimage.gaussian_filter(obj_sky.flux, smooth_sig)
#    else:
        #tmp = ndimage.gaussian_filter(arx_sky.flux, smooth_sig)
#        arx_skyspec = arx_skyspec.gauss_smooth(smooth_sig_pix*2*np.sqrt(2*np.log(2)))
        #arx_sky.flux = ndimage.gaussian_filter(arx_sky.flux, smooth_sig)

    # Define wavelengths of overlapping spectra
    keep_idx = np.where((obj_skyspec.wavelength.value>=min_wave) &
                         (obj_skyspec.wavelength.value<=max_wave))[0]
    #keep_wave = [i for i in obj_sky.wavelength.value if i>=min_wave if i<=max_wave]

    #Rebin both spectra onto overlapped wavelength range
    if len(keep_idx) <= 50:
        msgs.error("Not enough overlap between sky spectra")
    else: #rebin onto object ALWAYS
        keep_wave = obj_skyspec.wavelength[keep_idx]
        arx_skyspec = arx_skyspec.rebin(keep_wave)
        obj_skyspec = obj_skyspec.rebin(keep_wave)
        # Trim edges (rebinning is junk there)
        arx_skyspec.data['flux'][0,:2] = 0.
        arx_skyspec.data['flux'][0,-2:] = 0.
        obj_skyspec.data['flux'][0,:2] = 0.
        obj_skyspec.data['flux'][0,-2:] = 0.

    #   Normalize spectra to unit average sky count
    norm = np.sum(obj_skyspec.flux.value)/obj_skyspec.npix
    obj_skyspec.flux = obj_skyspec.flux / norm
    norm2 = np.sum(arx_skyspec.flux.value)/arx_skyspec.npix
    arx_skyspec.flux = arx_skyspec.flux / norm2
    if (norm < 0.):
        msgs.warn("Bad normalization of object in flexure algorithm")
        msgs.warn("Will try the median")
        norm = np.median(obj_skyspec.flux.value)
        if (norm < 0.):
            msgs.error("Improper sky spectrum for flexure.  Is it too faint??")
    if (norm2 < 0.):
        msgs.error("Bad normalization of archive in flexure. You are probably using wavelengths well beyond the archive.")

    #deal with bad pixels
    msgs.work("Need to mask bad pixels")

    #deal with underlying continuum
    msgs.work("Consider taking median first [5 pixel]")
    everyn = obj_skyspec.npix // 20
    mask, ct = arutils.robust_polyfit(obj_skyspec.wavelength.value, obj_skyspec.flux.value, 3, function='bspline',
                                  sigma=3., everyn=everyn)
    obj_sky_cont = arutils.func_val(ct, obj_skyspec.wavelength.value, 'bspline')
    obj_sky_flux = obj_skyspec.flux.value - obj_sky_cont
    mask, ct_arx = arutils.robust_polyfit(arx_skyspec.wavelength.value, arx_skyspec.flux.value, 3, function='bspline',
                                      sigma=3., everyn=everyn)
    arx_sky_cont = arutils.func_val(ct_arx, arx_skyspec.wavelength.value, 'bspline')
    arx_sky_flux = arx_skyspec.flux.value - arx_sky_cont

    # Consider shaprness filtering (e.g. LowRedux)
    msgs.work("Consider taking median first [5 pixel]")

    #Cross correlation of spectra
    #corr = np.correlate(arx_skyspec.flux, obj_skyspec.flux, "same")
    corr = np.correlate(arx_sky_flux, obj_sky_flux, "same")

    #Create array around the max of the correlation function for fitting for subpixel max
    # Restrict to pixels within maxshift of zero lag
    lag0 = corr.size//2
    mxshft = settings.argflag['reduce']['flexure']['maxshift']
    max_corr = np.argmax(corr[lag0-mxshft:lag0+mxshft]) + lag0-mxshft
    subpix_grid = np.linspace(max_corr-3., max_corr+3., 7.)

    #Fit a 2-degree polynomial to peak of correlation function
    fit = arutils.func_fit(subpix_grid, corr[subpix_grid.astype(np.int)], 'polynomial', 2)
    max_fit = -0.5*fit[1]/fit[2]

    #Calculate and apply shift in wavelength
    shift = float(max_fit)-lag0
    msgs.info("Flexure correction of {:g} pixels".format(shift))
    #model = (fit[2]*(subpix_grid**2.))+(fit[1]*subpix_grid)+fit[0]

    if msgs._debug['flexure']:
        debugger.plot1d(arx_skyspec.wavelength, arx_sky_flux, xtwo=np.roll(obj_skyspec.wavelength,int(-1*shift)), ytwo=obj_sky_flux)
        #debugger.xplot(arx_sky.wavelength, arx_sky.flux, xtwo=np.roll(obj_sky.wavelength.value,9), ytwo=obj_sky.flux*100)
        debugger.set_trace()

    flex_dict = dict(polyfit=fit, shift=shift, subpix=subpix_grid,
                     corr=corr[subpix_grid.astype(np.int)],
                     sky_spec=obj_skyspec,
                     arx_spec=arx_skyspec,
                     corr_cen=corr.size/2, smooth=smooth_sig_pix)
    # Return
    return flex_dict

