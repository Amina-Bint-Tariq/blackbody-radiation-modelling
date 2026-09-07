import numpy as np

from src.fitting import (
    fit_blackbody_spectrum,
    scaled_planck_radiance,
)


def test_fit_recovers_known_temperature():
    rng = np.random.default_rng(42)

    frequency = np.logspace(13.5, 15.5, 300)

    true_temperature = 6000.0
    true_scale = 2.5

    clean_spectrum = scaled_planck_radiance(
        frequency,
        true_scale,
        true_temperature,
    )

    uncertainty = 0.02 * np.max(clean_spectrum)

    noisy_spectrum = (
        clean_spectrum
        + rng.normal(
            0,
            uncertainty,
            size=frequency.size,
        )
    )

    uncertainties = np.full(
        frequency.shape,
        uncertainty,
    )

    result = fit_blackbody_spectrum(
        frequency,
        noisy_spectrum,
        uncertainties,
    )

    assert np.isclose(
        result["temperature"],
        true_temperature,
        rtol=0.05,
    )
