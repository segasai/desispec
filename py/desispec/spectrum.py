"""
desispec.spectrum
==============

Lightweight wrapper class for spectra, to be returned by io.read_frame
"""
from __future__ import absolute_import, division

import numpy as np

from desiutil.log import get_logger

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


