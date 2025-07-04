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


## 📊 Input Requirements

### 🧾 1. Prior Definitions

Before running the MCMC, the user is required to provide:



### 🔬 Physical and Model Parameters

| Parameter     | Description                                                                 |
|---------------|-----------------------------------------------------------------------------|
| `z`           | **Redshift** of the AGN.             |
| `theta`       | **Viewing angle** of the AGN in degrees. |
| `M`           | Initial guess and prior range for the **black hole mass** (in solar masses). |
| `Mdot`        | Initial guess and prior range for the **mass accretion rate** (in solar masses per year). |
| `logf`           | Initial guess and prior range for the **error term** used in the log-likelihood function. Represents fractional noise or model uncertainty. |

**Important to note that for the maximum Accretion Rate we use the physical constraint of the Eddington Limit**

### ⚙️ MCMC Configuration Parameters

| Parameter        | Description                                                                 |
|------------------|-----------------------------------------------------------------------------|
| `overlay_number` | Number of model spectra to overlay on the data plot from posterior samples. |
| `nwalkers`       | Number of walkers (independent chains) for the MCMC sampler. More walkers improve exploration but increase compute time. |
| `nsteps`         | Total number of steps each walker should take. Larger values result in better convergence. |
| `nburn`          | Number of initial steps to discard. These are excluded from the final posterior analysis. |

All of the parameters are entered **interactively** through prompts when the script runs!


### 🧪 Sample Input Values and Suggested Priors

Below is an example of the input values for a blazar:

| Parameter        | Example Value | Unit              | Suggested Prior Limits            | 
|------------------|---------------|-------------------|------------------------------------|
| `z`              | 3.41           | —                 | no prior                         | 
| `theta`          | 3             | degrees           | no prior                          | 
| `M`              | 2.5e9         | M☉                | 5e8 to 2e10                       | 
| `Mdot`           | 0.7           | M☉/year           | 1e-1 to Mdotedd                   | 
| `logf`           | -3            | —                  | -6 to -1                         | 
| `overlay_number` | 100           | —                 | no prior                          | 
| `nwalkers`       | 30            | —                 | no prior                          | 
| `nsteps`         | 2000          | —                 | no prior                          |
| `nburn`          | 900           | —                 | no prior                          | 

> These values are just examples. Adjust according to the specific AGN you are analyzing and your desired sampling precision.





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

## 📥 4.Required Python Packages

Before running the notebook, make sure the following Python packages are installed:

```
pip install numpy pandas scipy matplotlib emcee corner astropy ipython
```

Alternatively, install all dependencies with:

```
pip install -r requirements.txt
```

## 🚀 5.How to Run the Code

There are **two ways** to run the AGN accretion disk MCMC fitting code.

First method is used to integrate the script in bigger workflows and have flexibility and the second one for running the code independently


### 1️⃣Import and Use the MCMCAGNFit  Class


**Install Python & Required Packages**

Ensure Python 3.x is installed. Then install all required dependencies by downloading or cloning the repository and running:

```
git clone https://github.com/GeorgiosKyrio/mcmc-agnfit.git
```

and then

```
cd mcmc-agnfit
```

and then 

```
pip install -r requirements.txt
```


**Integration in the python script**

```
fit = MCMCAGNFit(
    csv_file=,
    z=,
    theta=,
    M_range=,
    Mdot_min=,
    logf_range=,
    initial_values=,
    overlay_number=,
    nwalkers=,
    nsteps=,
    nburn=
)

fit.run_sampler()   #Runs the MCMC algorithm and stores the samples.
fit.plot_corner()   #Creates and saves the corner plot of posterior distributions.
fit.plot_overlay()  #Plots observed data with multiple sampled model spectra overlaid
fit.summarize_posteriors()  #Calculates and saves the median and ±1σ intervals for each parameter.


```

**Again we have the Mdotedd as the physical constraint(maximum) for the accretion rate**


Write for example the following in your python script:

```
from mcmc_agnfit import MCMCAGNFit

fit = MCMCAGNFit(
    csv_file="data.csv",
    z=3.41,
    theta=3,
    M_range=[5e8, 2e10],
    Mdot_min=1e-1,
    logf_range=[-6, -1],
    initial_values=[2.5e9, 20, -3],
    overlay_number=100,
    nwalkers=30,
    nsteps=2000,
    nburn=900
)

fit.run_sampler()
fit.plot_corner()
fit.plot_overlay()
fit.summarize_posteriors()
```


> These values are just examples. Adjust according to the specific AGN you are analyzing and your desired sampling precision and **make sure that your data.csv file is in the same directory as your script and is in the appropriate form(see example)**.

### 2️⃣ Run the Script Directly (`mcmc-code.py`)

**Install Python & Required Packages**

Ensure Python 3.x is installed. Then install all required dependencies by downloading or cloning the repository and running:

```
git clone https://github.com/GeorgiosKyrio/mcmc-agnfit.git
```

and then

```
cd mcmc-agnfit
```

and then 

```
pip install -r requirements.txt
```

**Prepare Your Data File**

Have an appropriate CSV file(see example) named `data.csv` and place it in the **same directory** as `mcmc-code.py`.

**Run the Script**

Run the script and then input the necceasary parameters with:

```
python mcmc-code.py
```

---

## ✴️ Run With Dust Extinction Model (`MCMCAGNFitExt`)

The file `mcmc_agnfitext.py` provides an extended version of the AGN fitting model that incorporates **dust extinction** using the **Gordon et al. (2023)** parameterization of interstellar attenuation (G23). This allows for the fitting of **E(B–V)** as an additional physical parameter in the MCMC process (with Rv=3.1) .

### 🆕 Additional Model Parameter

| Parameter   | Description                             | Unit | Example Value | Typical Range |
|-------------|-----------------------------------------|------|----------------|----------------|
| `ebv_range` | Prior range for the color excess E(B–V) | mag  | (0.0, 0.5)      | (0.0 to 1.0)   |

This model uses the **G23 extinction law** implemented via `dust_extinction.parameter_averages.G23` to apply extinction at each frequency in the observer's frame.

---

### 📄 Data Format (Unchanged)

This extended version still expects input data in CSV format with the following columns:

```csv
frequency,freq_error,flux_freq,flux_freq_error
1.0e14,1.0e12,2.5e-13,1.0e-14

```

### ⚙️ Example: Using MCMCAGNFitExt with Extinction

```
from mcmc_agnfitext import MCMCAGNFitExt

fit = MCMCAGNFitExt(
    csv_file="data.csv",
    z=3.41,
    theta_deg=3,
    M_range=[5e8, 2e10],
    Mdot_min=1e-1,
    logf_range=[-6, -1],
    initial_values=[2.5e9, 20, -3, 0.05],  # includes E(B–V)
    overlay_number=100,
    nwalkers=30,
    nsteps=2000,
    nburn=900,
    ebv_range=(0.0, 0.5)  # user-defined prior on E(B–V)
)

fit.run_sampler()
fit.plot_corner()
fit.plot_overlay()
fit.summarize_posteriors()
```
---

## 📜 References

- Shakura, N. I., & Sunyaev, R. A. (1973). *Black holes in binary systems. Observational appearance.* Astronomy and Astrophysics, **24**, 337–355.
- Foreman-Mackey, D., Hogg, D. W., Lang, D., & Goodman, J. (2013). *emcee: The MCMC Hammer.* Publications of the Astronomical Society of the Pacific, **125**(925), 306–312. DOI: [10.1086/670067](https://doi.org/10.1086/670067)
- Gordon, K. D., Clayton, G. C., Misselt, K. A., & Wolf, A. (2023). A New Parameterization of the Interstellar Extinction Curve from the Ultraviolet to the Near-infrared. Astrophysical Journal. DOI: 10.3847/1538-4357/acf3e3

---
