# mcmc-agnfit
Spectral fitting of AGN data using MCMC sampling for astrophysical parameter estimation.
# MCMC-AGNFit

**mcmc-agnfit** is a Python tool for fitting Active Galactic Nucleus (AGN) spectra using the Shakura–Sunyaev accretion disk model and Markov Chain Monte Carlo (MCMC) sampling via the `emcee` package.

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

- **Redshift** `z` of the AGN under investigation
- **Viewing angle** `theta` of the AGN in degrees
- **Initial guesses** and **prior ranges** for the following parameters:
  - `log_M` – Logarithmic black hole mass (base 10), in solar masses
  - `log_Mdot` – Logarithmic accretion rate, in solar masses per year
  - `log_f` – Logarithmic error used in the log_likelihood function

These priors and initials are input interactively or via a cell in the notebook.
- **overlay_number** `n` that gives the number of spectra that are going to be shown in the overlay graph

### 📄 2. CSV Data Format

Users must supply a CSV file named `data.csv` with the following **required columns**:

| Column               | Description                          |
|----------------------|--------------------------------------|
| `frequency`          | Frequency in Hz                      |
| `freq_error`         | Frequency error                      |
| `flux_freq`          | νFν value in CGS units               |
| `flux_freq_error`    | Uncertainty in νFν                   |

Example CSV:

```csv
frequency,freq_error,flux_freq,flux_freq_error
1.0e14,1.0e12,2.5e-13,1.0e-14
```

## 📦 3. Output Files

After running the notebook, the following output files are generated:

| File              | Description                                                                 |
|-------------------|-----------------------------------------------------------------------------|
| `output.h5`       | HDF5 file containing the raw MCMC chains for the sampled parameters         |
| `corner_plot.png` | Corner plot showing the posterior distributions of the model parameters     |
| `overlay`         | Overlay of the accretion disk spectra for different posterior parameters    |
| *Terminal Output* | A printed summary of each parameter's median and uncertainty range          |

### 📈 Posterior Summary

The MCMC chains produce the posterior distributions for each model parameter. From these distributions, we extract:

- The **median value** (50th percentile), $\tilde{\theta} = p_{50}$,
- The **lower bound** of the 1σ confidence interval, $\Delta_{-} = p_{50} - p_{16}$,
- The **upper bound** of the 1σ confidence interval, $\Delta_{+} = p_{84} - p_{50}$,

where $p_{16}$, $p_{50}$, and $p_{84}$ are the 16th, 50th, and 84th percentiles of the sampled posterior distribution.

## 📥 Required Python Packages

Before running the notebook, make sure the following Python packages are installed:

```
pip install numpy pandas scipy matplotlib emcee corner astropy h5py
```

Alternatively, install all dependencies with:

```
pip install -r requirements.txt
```

