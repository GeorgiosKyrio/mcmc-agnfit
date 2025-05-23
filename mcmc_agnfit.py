import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import pandas as pd
import emcee
from astropy.io import ascii
import corner
from multiprocessing import Pool
import os

# Constants in CGS units
c = 3e10
h = 6.626e-27
k_B = 1.38e-16
sigma = 5.67e-5
G = 6.674e-8
M_sun = 1.989e33
L_Edd_factor = 1.26e38
H0 = 2.195e-18
Omega_m = 0.308
Omega_Lambda = 0.692


class MCMCAGNFit:
    def __init__(self, csv_file, z, theta, M_range, Mdot_min, logf_range,
                 initial_values, overlay_number, nwalkers, nsteps, nburn):
        self.csv_file = csv_file
        self.z = z
        self.theta = np.deg2rad(theta)
        self.M_range = M_range
        self.Mdot_min = Mdot_min
        self.logf_range = logf_range
        self.initial_values = initial_values
        self.overlay_number = overlay_number
        self.nwalkers = nwalkers
        self.nsteps = nsteps
        self.nburn = nburn
        self.frequencies = np.logspace(10, 18, 100)
        self.load_data()

    def load_data(self):
        data = pd.read_csv(self.csv_file)
        self.freq_log = np.log10(data.iloc[:, 0].values)
        self.flux_log = np.log10(data.iloc[:, 1].values)
        flux_err = 0.1 * data.iloc[:, 1].values
        self.flux_log_err = np.log10(data.iloc[:, 1].values + flux_err) - np.log10(data.iloc[:, 1].values)

    def luminosity_distance(self, z):
        def integrand(z_prime):
            return 1.0 / np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_Lambda)
        integral, _ = quad(integrand, 0, z)
        return c * (1 + z) * integral / H0

    def disk_temperature(self, r, M_bh, M_dot):
        R_s = 2 * G * M_bh / c**2
        R_in = 3 * R_s
        if r < R_in:
            return 0.0
        factor = 3 * M_dot * R_s * c**2 / (16 * np.pi * sigma * r**3)
        T = (factor * (1 - np.sqrt(R_in / r)))**0.25
        return np.clip(T, 1e1, 1e12)

    def planck_function(self, nu, T):
        if T == 0.0:
            return 0.0
        exponent = h * nu / (k_B * T)
        if exponent > 700:
            return 0.0
        return (2 * h * nu**3 / c**2) / (np.exp(exponent) - 1)

    def flux_density(self, nu, M_bh, M_dot):
        R_s = 2 * G * M_bh / c**2
        R_in = 3 * R_s
        R_out = 1000 * R_s
        def integrand(r):
            return 2 * np.pi * r * np.cos(self.theta) * self.planck_function(nu, self.disk_temperature(r, M_bh, M_dot))
        r_vals = np.logspace(np.log10(R_in), np.log10(R_out), 100)
        integrand_vals = [integrand(r) for r in r_vals]
        integral = np.trapezoid(np.array(integrand_vals)*r_vals, np.log(r_vals))
        d_L = self.luminosity_distance(self.z)
        return integral / d_L**2

    def simulate_spectrum(self, freqs, M_bh, M_dot):
        return np.array([self.flux_density(nu, M_bh, M_dot) for nu in freqs * (1 + self.z)])

    def log_prior(self, params):
        M, Mdot, logf = params
        
        M_linear=10**M
        Mdot_linear=10**Mdot
        
        Mdot_Edd = L_Edd_factor * (M_linear / M_sun) / (0.1 * c**2)
        if np.log10(self.M_range[0])+np.log10(M_sun) < M < np.log10(self.M_range[1])+np.log10(M_sun) and \
           np.log10(self.Mdot_min)+np.log10(M_sun/(365*24*3600)) < Mdot < np.log10(Mdot_Edd) and \
           self.logf_range[0] < logf < self.logf_range[1]:
            return 0.0
        return -np.inf

    def log_likelihood(self, params):
        M, Mdot, logf = params

        M_linear=10**M
        Mdot_linear=10**Mdot
        
        model_fluxes = self.simulate_spectrum(10**self.freq_log, M_linear, Mdot_linear) * (10**self.freq_log) * (1 + self.z)
        model_log = np.log10(model_fluxes)
        sigma2 = self.flux_log_err**2 + np.exp(2 * logf)
        return -0.5 * np.sum((self.flux_log - model_log)**2 / sigma2 + np.log(sigma2))

    def log_posterior(self, params):
        lp = self.log_prior(params)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(params)

    def run_sampler(self):
        ndim = 3
        initial = [np.log10(self.initial_values[0] * M_sun),
                   np.log10(self.initial_values[1] * M_sun / (365*24*3600)),
                   self.initial_values[2]]
        pos = initial + np.abs(2e-1 * np.random.randn(self.nwalkers, ndim))

        out_dir = "./"  # or any directory you want
        filename = os.path.join(out_dir, "output.h5")

        backend = emcee.backends.HDFBackend("output.h5")
        backend.reset(self.nwalkers, ndim)

        with Pool() as pool:
            sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.log_posterior, pool=pool, backend=backend)
            sampler.run_mcmc(pos, self.nsteps, progress=True)
        
        
        
        self.samples = sampler.get_chain(discard=self.nburn, thin=15, flat=True)

    def plot_corner(self):
        M_bh = 10**self.samples[:, 0] / M_sun
        Mdot = 10**self.samples[:, 1] * 365*24*3600 / M_sun
        logf = self.samples[:, 2]

        fig = corner.corner(np.vstack([M_bh, Mdot, logf]).T,
                            labels=["$M_{bh}$ ($M_\odot$)", "$\\dot{M}$ ($M_\odot$/yr)", "log$_{10}(f)$"],
                            truths=self.initial_values)
        fig.savefig("corner_plot.png")

    def plot_overlay(self):
        freqs = self.frequencies
        inds = np.random.randint(len(self.samples), size=self.overlay_number)

        plt.figure(figsize=(8, 6))
        for ind in inds:
            M = 10**self.samples[ind, 0] / M_sun
            Mdot = 10**self.samples[ind, 1] * 365*24*3600 / M_sun
            M_cgs = M * M_sun
            Mdot_cgs = Mdot * M_sun / (365*24*3600)
            spectrum = self.simulate_spectrum(freqs, M_cgs, Mdot_cgs)
            plt.loglog(freqs, spectrum * freqs * (1 + self.z), color='orange', alpha=0.3)

        plt.loglog(10**self.freq_log, 10**self.flux_log, 'bo', label="Observed Data")
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Flux Density * Frequency [erg/s/cm²]")
        plt.ylim(1e-15,1e-11)
        plt.xlim(1e12,1e16)
        plt.title("Overlay of Model Spectra")
        plt.grid(True, which="both", ls="--")
        plt.savefig("model_overlay.png")

    def summarize_posteriors(self):
        flat_samples = self.samples.copy()
        flat_samples[:, 0:2] = 10**flat_samples[:, 0:2]
        flat_samples[:, 0] /= M_sun
        flat_samples[:, 1] *= 365 * 24 * 3600 / M_sun

        summary = {}
        for i, name in enumerate(["M_bh", "Mdot", "logf"]):
            p16, p50, p84 = np.percentile(flat_samples[:, i], [16, 50, 84])
            summary[name] = {
                "median": p50,
                "-1sigma": p50 - p16,
                "+1sigma": p84 - p50
            }

        df = pd.DataFrame(summary).T
        df.to_csv("fit_results.csv")
        print(df)
        return df
