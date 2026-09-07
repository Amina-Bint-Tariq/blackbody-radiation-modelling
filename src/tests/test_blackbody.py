import numpy as np

from scipy.constants import Stefan_Boltzmann as sigma

from src.blackbody import (
    peak_frequency,
    planck_radiance,
    total_luminosity_numerical,
)


def test_planck_radiance_is_positive():
    frequency = 3e14
    temperature = 5000

    result = planck_radiance(frequency, temperature)

    assert result > 0


def test_hotter_blackbody_peaks_at_higher_frequency():
    cool_peak = peak_frequency(4000)
    hot_peak = peak_frequency(8000)

    assert hot_peak > cool_peak


def test_peak_frequency_scales_with_temperature():
    peak_4000 = peak_frequency(4000)
    peak_8000 = peak_frequency(8000)

    ratio = peak_8000 / peak_4000

    assert np.isclose(ratio, 2.0, rtol=0.01)


def test_numerical_luminosity_matches_stefan_boltzmann():
    radius = 1.0
    temperature = 5000
    emissivity = 1.0

    numerical = total_luminosity_numerical(
        radius,
        temperature,
        emissivity,
    )

    expected = (
        4
        * np.pi
        * radius**2
        * emissivity
        * sigma
        * temperature**4
    )

    assert np.isclose(
        numerical,
        expected,
        rtol=0.02,
    )


def test_negative_temperature_raises_error():
    try:
        planck_radiance(3e14, -5000)

    except ValueError:
        pass

    else:
        raise AssertionError(
            "Negative temperature should raise ValueError."
        )
