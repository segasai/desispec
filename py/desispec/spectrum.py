"""
desispec.spectrum
==============

Lightweight wrapper class for spectra, to be returned by io.read_frame
"""
from __future__ import absolute_import, division

import numpy as np

from desispec.interpolation import resample_flux

class Spectrum(object):
    def __init__(self, wave, flux, ivar, mask=None, R=None):
        """Lightweight wrapper of a single spectrum

        Args:
            wave (1D ndarray): wavelength in Angstroms
            flux (1D ndarray): flux (photons or ergs/s/cm^2/A)
            ivar (1D ndarray): inverse variance of flux
            R : Resolution object

        All args become attributes.  This is syntactic sugar.
        """
        self.wave = wave
        self.flux = flux
        self.ivar = ivar
        if mask is None:
            self.mask = np.zeros(self.flux.shape, dtype=int)
        else:
            self.mask = mask

        self.R = R

    @property
    def npix(self):
        return len(self.wave)

    @property
    def sig(self):
        # Sigma from ivar
        sig = np.zeros_like(self.ivar)
        gd_ivar = self.ivar > 0.
        sig[gd_ivar] = 1./np.sqrt(self.ivar[gd_ivar])
        # Return
        return sig

    def plot(self):
        from matplotlib import pyplot as plt
        # Simple plot
        plt.clf()
        ax = plt.gca()
        ax.plot(self.wave, self.flux)
        ax.plot(self.wave, self.sig, color='red')
        plt.show()

    def rebin(self, new_wave):
        """ Rebin the current spectrum onto a new wavelength array
        Uses desispec.interpolation.resample_flux

        Args:
            new_wave:  ndarray
              New wavelength values

        Returns:
        new_spec : Spectrum

        """
        new_flux, new_ivar = resample_flux(new_wave, self.wave, self.flux, ivar=self.ivar)
        # Return
        return Spectrum(new_wave, new_flux, new_flux)

