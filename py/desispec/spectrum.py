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


