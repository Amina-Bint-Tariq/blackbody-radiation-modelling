"""
Parameter-estimation utilities for black-body spectra.

This module contains functions for fitting a scaled Planck spectrum
to spectral measurements and extracting uncertainties from the
resulting covariance matrix.
"""

import numpy as np

from scipy.constants import Boltzmann as k_B
from scipy.constants import Planck as h
from scipy.optimize import curve_fit

from .blackbody import planck_radiance


WIEN_FREQUENCY_CONSTANT = 2.821439


def scaled_planck_radiance(frequency, scale, temperature):
    """
    Calculate a scaled black-body spectrum.

    The scale parameter absorbs geometric factors such as source size,
    distance and instrumental calibration.

    Parameters
    ----------
    frequency : float or array-like
        Frequency in Hz.
    scale : float
        Multiplicative scale factor.
    temperature : float
        Temperature in K.

    Returns
    -------
    float or numpy.ndarray
        Scaled spectral radiance.
    """

    return scale * planck_radiance(frequency, temperature)


def initial_temperature_guess(frequency, intensity):
    """
    Estimate temperature from the observed peak frequency.

    Parameters
    ----------
    frequency : array-like
        Frequencies in Hz.
    intensity : array-like
        Observed spectral values.

    Returns
    -------
    float
        Initial temperature estimate in K.
    """

    frequency = np.asarray(frequency, dtype=float)
    intensity = np.asarray(intensity, dtype=float)

    if frequency.shape != intensity.shape:
        raise ValueError(
            "Frequency and intensity arrays must have matching shapes."
        )

    peak_index = np.argmax(intensity)
    peak_frequency = frequency[peak_index]

    temperature_guess = (
        h * peak_frequency
        / (WIEN_FREQUENCY_CONSTANT * k_B)
    )

    return temperature_guess


def fit_blackbody_spectrum(
    frequency,
    intensity,
    uncertainty=None,
):
    """
    Fit a scaled Planck spectrum to observations.

    Parameters
    ----------
    frequency : array-like
        Observed frequencies in Hz.
    intensity : array-like
        Observed spectral values.
    uncertainty : array-like, optional
        One-sigma uncertainty associated with each spectral value.

    Returns
    -------
    dict
        Dictionary containing the fitted temperature, scale factor,
        their uncertainties and the full covariance matrix.
    """

    frequency = np.asarray(frequency, dtype=float)
    intensity = np.asarray(intensity, dtype=float)

    if frequency.shape != intensity.shape:
        raise ValueError(
            "Frequency and intensity arrays must have matching shapes."
        )

    if uncertainty is not None:
        uncertainty = np.asarray(uncertainty, dtype=float)

        if uncertainty.shape != intensity.shape:
            raise ValueError(
                "Uncertainty must have the same shape as intensity."
            )

        if np.any(uncertainty <= 0):
            raise ValueError(
                "All uncertainty values must be greater than zero."
            )

    temperature_guess = initial_temperature_guess(
        frequency,
        intensity,
    )

    peak_index = np.argmax(intensity)

    unit_model_at_peak = planck_radiance(
        frequency[peak_index],
        temperature_guess,
    )

    scale_guess = (
        intensity[peak_index]
        / unit_model_at_peak
    )

    popt, pcov = curve_fit(
        scaled_planck_radiance,
        frequency,
        intensity,
        p0=(scale_guess, temperature_guess),
        sigma=uncertainty,
        absolute_sigma=uncertainty is not None,
        maxfev=20_000,
    )

    scale, temperature = popt

    parameter_uncertainties = np.sqrt(np.diag(pcov))

    scale_uncertainty = parameter_uncertainties[0]
    temperature_uncertainty = parameter_uncertainties[1]

    return {
        "temperature": temperature,
        "temperature_uncertainty": temperature_uncertainty,
        "scale": scale,
        "scale_uncertainty": scale_uncertainty,
        "covariance": pcov,
    }
