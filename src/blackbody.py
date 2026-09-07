"""
Physical models and numerical utilities for black-body radiation.

This module implements Planck's law in frequency space and provides
functions for calculating the spectral luminosity, total luminosity,
and peak-emission frequency of an idealised spherical black body.
"""

import numpy as np

from scipy.constants import Boltzmann as k_B
from scipy.constants import Planck as h
from scipy.constants import c
from scipy.integrate import quad
from scipy.optimize import minimize_scalar


def planck_radiance(frequency, temperature):
    """
    Calculate black-body spectral radiance using Planck's law.

    Parameters
    ----------
    frequency : float or array-like
        Frequency in Hz.
    temperature : float
        Black-body temperature in K.

    Returns
    -------
    float or numpy.ndarray
        Spectral radiance B_nu in W m^-2 sr^-1 Hz^-1.

    Notes
    -----
    Planck's law in frequency space is

        B_nu(T) = (2 h nu^3 / c^2) /
                  (exp(h nu / (k_B T)) - 1)
    """

    frequency = np.asarray(frequency, dtype=float)

    if temperature <= 0:
        raise ValueError("Temperature must be greater than zero.")

    if np.any(frequency < 0):
        raise ValueError("Frequency values must be non-negative.")

    exponent = h * frequency / (k_B * temperature)

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        radiance = (
            2 * h * frequency**3 / c**2
        ) / np.expm1(exponent)

    return np.nan_to_num(radiance)


def spectral_luminosity(frequency, radius, temperature, emissivity=1.0):
    """
    Calculate the spectral luminosity of a spherical black body.

    Parameters
    ----------
    frequency : float or array-like
        Frequency in Hz.
    radius : float
        Radius of the emitting sphere in metres.
    temperature : float
        Temperature in K.
    emissivity : float, optional
        Dimensionless emissivity between 0 and 1.
        Default is 1.0.

    Returns
    -------
    float or numpy.ndarray
        Spectral luminosity in W Hz^-1.
    """

    if radius <= 0:
        raise ValueError("Radius must be greater than zero.")

    if not 0 <= emissivity <= 1:
        raise ValueError("Emissivity must lie between 0 and 1.")

    surface_area = 4 * np.pi * radius**2

    # Surface flux density is pi * B_nu for isotropic emission.
    spectral_flux = np.pi * planck_radiance(frequency, temperature)

    return surface_area * emissivity * spectral_flux


def total_luminosity_numerical(
    radius,
    temperature,
    emissivity=1.0,
    upper_frequency=1e16,
):
    """
    Numerically integrate the spectral luminosity.

    Parameters
    ----------
    radius : float
        Radius of the black body in metres.
    temperature : float
        Temperature in K.
    emissivity : float, optional
        Dimensionless emissivity between 0 and 1.
    upper_frequency : float, optional
        Upper integration limit in Hz.

    Returns
    -------
    float
        Total emitted luminosity in W.
    """

    if upper_frequency <= 0:
        raise ValueError("Upper frequency must be greater than zero.")

    luminosity, _ = quad(
        lambda nu: spectral_luminosity(
            nu,
            radius,
            temperature,
            emissivity,
        ),
        0,
        upper_frequency,
        limit=500,
    )

    return luminosity


def peak_frequency(temperature):
    """
    Numerically determine the frequency at which Planck radiance peaks.

    Parameters
    ----------
    temperature : float
        Temperature in K.

    Returns
    -------
    float
        Peak frequency in Hz.
    """

    if temperature <= 0:
        raise ValueError("Temperature must be greater than zero.")

    # Wien-like estimate in frequency space gives a good search centre.
    initial_peak = 2.821439 * k_B * temperature / h

    lower = initial_peak / 10
    upper = initial_peak * 10

    result = minimize_scalar(
        lambda nu: -planck_radiance(nu, temperature),
        bounds=(lower, upper),
        method="bounded",
    )

    if not result.success:
        raise RuntimeError("Peak-frequency optimisation failed.")

    return result.x
