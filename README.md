# mcmc-agnfit
Spectral fitting of AGN data using MCMC sampling for astrophysical parameter estimation.
# MCMC-AGNFit

**MCMC-AGNFit** is a Python tool for fitting Active Galactic Nucleus (AGN) spectra using the Shakura–Sunyaev accretion disk model and Markov Chain Monte Carlo (MCMC) sampling via the `emcee` package.

---

## 🔭 Scientific Background

The code models AGN continuum emission using a standard **Shakura–Sunyaev thin accretion disk** model, a widely used description of disk emission around supermassive black holes. The spectral energy distribution is computed based on physical parameters:

- **Black hole mass** (M)
- **Mass accretion rate** (ṁ)
- **Viewing angle (cos θ)** 

The disk spectrum is compared to observed photometric data in frequency units, using Bayesian inference to estimate posterior distributions for these physical parameters.

---

## 📜 References

- Shakura, N. I., & Sunyaev, R. A. (1973). *Black holes in binary systems. Observational appearance.* Astronomy and Astrophysics, **24**, 337–355.
- Foreman-Mackey, D., Hogg, D. W., Lang, D., & Goodman, J. (2013). *emcee: The MCMC Hammer.* Publications of the Astronomical Society of the Pacific, **125**(925), 306–312. DOI: [10.1086/670067](https://doi.org/10.1086/670067)

---

## 📊 Input Requirements

### 🧾 1. Prior Definitions

Before running the MCMC, the user is required to provide:

-redshift `z` of the AGN under investigation
-viewing angle theta of the AGN in degrees

- **Initial guesses** and **prior ranges** for the following parameters:
  - `log_M` – Logarithmic black hole mass (base 10), in solar masses
  - `log_Mdot` – Logarithmic accretion rate, in solar masses per year
  - `log_f` – Logarithmic error used in the log_likelihood function

These priors and initials are input interactively or via a cell in the notebook.

### 📄 2. CSV Data Format

Users must supply a CSV file with the following **required columns**:

| Column               | Description                          |
|----------------------|--------------------------------------|
| `frequency`          | Frequency in Hz                      |
| `freq_error`         | Frequency error (optional)           |
| `flux_freq`          | νFν value in CGS units               |
| `flux_freq_error`    | Uncertainty in νFν                   |

Example CSV:

```csv
frequency,freq_error,flux_freq,flux_freq_error
1.0e14,1.0e12,2.5e-13,1.0e-14
...
