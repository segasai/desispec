"""
desispec.flexure
==================

Functions to measure 'flexure' in the instrument using
the night sky lines

A fair bit of this has been ported/grabbed from PYPIT
"""
from __future__ import print_function, absolute_import, division

import numpy as np
from pkg_resources import resource_exists, resource_filename

from astropy.io import fits

from desiutil.log import get_logger
from desiutil.funcfits import iter_fit, func_fit, func_val
from desispec.spectrum import Spectrum
from desispec.bootcalib import find_arc_lines


def load_paranal():
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


def flex_shift(channel, obj_skyspec, mxshft=10., debug=False):
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
    mxshft : float
      Maximum shift to allow for in the flexure analysis
    dwv_arx : float
      Dispersion of Paranal spectrum
      R_Paranal = 31300

    Returns
    -------
    flex_dict
    """
    from astropy.convolution import convolve, Gaussian1DKernel

    # Load archive spectrum
    arx_skyspec = load_paranal()

    # Logging
    log=get_logger()
    log.info("starting flexure shift calculation")
    # DESI resolution (sigma)
    wdisp_dict = dict(b=0.628, r=0.578, z=0.731)  # Ang
    wvmed_dict = dict(b=4765., r=6685., z=8635.) # Ang
    # Paranal (FWHM)
    R_Paranal = 31300.
    wdisp_paranal = wvmed_dict[channel] / R_Paranal

    # Find the brightest emission lines
    arx_cent, arx_amp = find_arc_lines(arx_skyspec.flux)
    obj_cent, obj_amp = find_arc_lines(obj_skyspec.flux)

    '''
    #Keep only 5 brightest amplitude lines (xxx_keep is array of indices within arx_w of the 5 brightest)
    arx_keep = np.argsort(arx_amp)[-5:]
    obj_keep = np.argsort(obj_amp)[-5:]

    '''

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

    # Calculate dispersion (Angstrom per pixel)
    arx_disp = np.append(arx_skyspec.wave[1]-arx_skyspec.wave[0],
                         arx_skyspec.wave[1:]-arx_skyspec.wave[:-1])

    # Determine sigma of gaussian for smoothing
    arx_med_sig2 = (wdisp_paranal/(2*np.sqrt(2*np.log(2))))**2
    obj_med_sig2 = wdisp_dict[channel]**2

    # Smooth
    smooth_sig = np.sqrt(obj_med_sig2-arx_med_sig2)  # Ang
    smooth_sig_pix = smooth_sig / np.median(arx_disp)

    # gaussian drops to 1/100 of maximum value at x =
    # sqrt(2*ln(100))*sigma, so number of pixels to include from
    # centre of gaussian is:
    const100 = 3.034854259             # sqrt(2*ln(100))
    n = np.ceil(const100 * smooth_sig_pix)
    x_size = int(2*n) + 1  # we want this to be odd integer
    arx_skyspec.flux = convolve(arx_skyspec.flux,
                                Gaussian1DKernel(smooth_sig_pix, x_size=x_size),
                                boundary='fill', fill_value=0., normalize_kernel=True)

    # Determine region of wavelength overlap
    min_wave = max(np.amin(arx_skyspec.wave), np.amin(obj_skyspec.wave))
    max_wave = min(np.amax(arx_skyspec.wave), np.amax(obj_skyspec.wave))

    # Define wavelengths of overlapping spectra
    keep_idx = np.where((obj_skyspec.wave >= min_wave) &
                         (obj_skyspec.wave <= max_wave))[0]

    # Rebin both spectra onto overlapped wavelength range
    # rebin onto object ALWAYS
    keep_wave = obj_skyspec.wave[keep_idx]
    arx_skyspec = arx_skyspec.rebin(keep_wave)
    obj_skyspec = obj_skyspec.rebin(keep_wave)
    # Trim edges (rebinning is junk there)
    arx_skyspec.flux[:2] = 0.
    arx_skyspec.flux[-2:] = 0.
    obj_skyspec.flux[:2] = 0.
    obj_skyspec.flux[-2:] = 0.

    #   Normalize spectra to unit average sky count
    norm = np.sum(obj_skyspec.flux)/obj_skyspec.npix
    obj_skyspec.flux = obj_skyspec.flux / norm
    norm2 = np.sum(arx_skyspec.flux)/arx_skyspec.npix
    arx_skyspec.flux = arx_skyspec.flux / norm2
    if (norm < 0.):
        log.warning("Bad normalization of object in flexure algorithm")
        log.warning("Will try the median")
        norm = np.median(obj_skyspec.flux.value)
        if (norm < 0.):
            log.error("Improper sky spectrum for flexure.  Is it too faint??")
    if (norm2 < 0.):
        log.error("Bad normalization of archive in flexure. You are probably using wavelengths well beyond the archive.")

    # Deal with underlying continuum
    everyn = obj_skyspec.npix // 20
    obj_fdict, _ = iter_fit(obj_skyspec.wave, obj_skyspec.flux, 'bspline', 3, sig_rej=3., everyn=everyn)
    obj_sky_cont = func_val(obj_skyspec.wave, obj_fdict)
    obj_sky_flux = obj_skyspec.flux - obj_sky_cont
    arx_fdict, mask = iter_fit(arx_skyspec.wave, arx_skyspec.flux, 'bspline', 3, sig_rej=3., everyn=everyn)
    arx_sky_cont = func_val(arx_skyspec.wave, arx_fdict)
    arx_sky_flux = arx_skyspec.flux - arx_sky_cont

    #Cross correlation of spectra
    corr = np.correlate(arx_sky_flux, obj_sky_flux, "same")

    #Create array around the max of the correlation function for fitting for subpixel max
    # Restrict to pixels within maxshift of zero lag
    lag0 = corr.size//2
    max_corr = np.argmax(corr[lag0-mxshft:lag0+mxshft]) + lag0-mxshft
    subpix_grid = np.linspace(max_corr-3., max_corr+3., 7.)

    #Fit a 2-degree polynomial to peak of correlation function
    fit = func_fit(subpix_grid, corr[subpix_grid.astype(np.int)], 'polynomial', 2)
    max_fit = -0.5*fit[1]/fit[2]

    #Calculate and apply shift in wavelength
    shift = float(max_fit)-lag0
    log.info("Flexure correction of {:g} pixels".format(shift))

    if debug:
        import pdb; pdb.set_trace()
        #debugger.plot1d(arx_skyspec.wavelength, arx_sky_flux, xtwo=np.roll(obj_skyspec.wavelength,int(-1*shift)), ytwo=obj_sky_flux)
        #debugger.set_trace()

    flex_dict = dict(polyfit=fit, shift=shift, subpix=subpix_grid,
                     corr=corr[subpix_grid.astype(np.int)],
                     sky_spec=obj_skyspec,
                     arx_spec=arx_skyspec,
                     corr_cen=corr.size/2, smooth=smooth_sig_pix)
    # Return
    return flex_dict

