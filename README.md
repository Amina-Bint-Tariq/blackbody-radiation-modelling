# Black-Body Radiation Modelling

A computational physics project exploring the numerical modelling and
parameter estimation of black-body radiation using Python.

## Project overview

This project uses numerical methods to investigate black-body radiation,
including:

- modelling spectra using Planck's law
- numerical integration of emitted radiation
- numerical determination of spectral peaks
- fitting spectral measurements to estimate temperature
- uncertainty estimation and residual analysis

## Tools

- Python
- NumPy
- SciPy
- Matplotlib
- Jupyter

## Status

Project currently being developed and expanded from previous university
scientific-programming work.
## Running the project locally

Clone the repository and create a virtual environment:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

Install the project dependencies:

```powershell
python -m pip install -r requirements.txt
```

Run the automated test suite:

```powershell
python -m pytest -q
```

The tests validate the numerical black-body calculations and confirm that the fitting routine can recover the temperature of a synthetic spectrum.