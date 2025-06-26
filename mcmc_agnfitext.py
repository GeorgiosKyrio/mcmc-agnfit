import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
import emcee
from dust_extinction.parameter_averages import G23
import astropy.units as u
from multiprocessing import Pool
import os

# Constants (CGS units unless noted)
c = 3e10            # speed of light [cm/s]
h = 6.626e-27       # Planck constant [erg·s]
k_B = 1.38e-16      # Boltzmann constant [erg/K]
sigma = 5.67e-5     # Stefan-Boltzmann [erg/cm^2/s/K^4]
G = 6.674e-8        # gravitational constant [cm^3/g/s^2]
M_sun = 1.989e33    # solar mass [g]
L_Edd_factor = 1.26e38  # Eddington luminosity factor [erg/s]
H0 = 2.195e-18      # Hubble constant [s^-1]
Omega_m = 0.308     # cosmological Omega_matter
Omega_Lambda = 0.692# cosmological Omega_lambda
Rv = 3.1            # total-to-selective extinction ratio

class MCMCAGNFitExt:
    def __init__(self, csv_file, z, theta_deg, M_range, Mdot_min, logf_range,
                 initial_values, overlay_number, nwalkers, nsteps, nburn,
                 ebv_range):
        """
        ...
        ebv_range: tuple (min, max) range for E(B-V)
        """
        self.csv_file = csv_file
        self.z = z
        self.theta = np.deg2rad(theta_deg)
        self.M_range = (np.log10(M_range[0]*M_sun), np.log10(M_range[1]*M_sun))
        self.Mdot_min = np.log10(Mdot_min * M_sun / (365*24*3600))
        self.logf_range = logf_range
        self.initial_values = initial_values
        self.overlay_number = overlay_number
        self.nwalkers = nwalkers
        self.nsteps = nsteps
        self.nburn = nburn
        self.ebv_range = ebv_range

        self.frequencies = np.logspace(13, 16, 100)
        self.load_data()

     def load_data(self):
        data = pd.read_csv(self.csv_file)
        self.freq_log = np.log10(data.iloc[:, 0].values)
        self.flux_log = np.log10(data.iloc[:, 1].values)
        flux_err = 0.1 * data.iloc[:, 1].values
        self.flux_log_err = np.log10(data.iloc[:, 1].values + flux_err) - np.log10(data.iloc[:, 1].values)

    def luminosity_distance(self, z):
        """Compute luminosity distance [cm] for redshift z (flat LambdaCDM)."""
        def integrand(zp):
            return 1.0 / np.sqrt(Omega_m * (1 + zp)**3 + Omega_Lambda)
        integral, _ = quad(integrand, 0, z)
        return c * (1 + z) * integral / H0

    def disk_temperature(self, r, M_bh, M_dot):
        """Temperature profile of a thin accretion disk."""
        R_s = 2 * G * M_bh / c**2
        R_in = 3 * R_s
        if r < R_in:
            return 0.0
        factor = 3 * M_dot * R_s * c**2 / (16 * np.pi * sigma * r**3)
        T = (factor * (1 - np.sqrt(R_in/r)))**0.25
        return np.clip(T, 1e1, 1e12)

    def planck_function(self, nu, T):
        """Blackbody Planck function at frequency nu and temperature T."""
        if T <= 0.0:
            return 0.0
        exponent = h * nu / (k_B * T)
        if exponent > 700:  # avoid overflow
            return 0.0
        return (2 * h * nu**3 / c**2) / (np.exp(exponent) - 1)

    def extinction_factor(self, nu, ebv):
        """
        Attenuation factor for frequency nu (Hz) given E(B-V) [mag].
        Uses the G23 extinction curve with Rv = 3.1:contentReference[oaicite:3]{index=3}:contentReference[oaicite:4]{index=4}.
        """
        # Convert frequency to wavelength in microns
        wavelength_um = (c / nu) * 1e4  # cm -> micron
        x_inv_micron = 1.0 / wavelength_um  # [1/micron]
        # G23 expects input in units of 1/micron:
        ext_model = G23(Rv=Rv)  # Gordon 2023 model
        # Evaluate A(λ)/A(V) using G23 model
        # Note: G23 can take astropy units or numeric array with units.
        A_over_Av = ext_model(x_inv_micron * u.Unit("1/um"))
        # A(V) = Rv * E(B-V), so A(λ) = A(V) * [A(λ)/A(V)]
        A_lambda = (Rv * ebv) * A_over_Av
        # Attenuation factor = 10^{-0.4 A(λ)}
        return 10.0 ** (-0.4 * A_lambda.value)  # .value to strip units

    def flux_density(self, nu, M_bh, M_dot, ebv):
        """Compute disk flux density F_nu at frequency nu, including extinction."""
        R_s = 2 * G * M_bh / c**2
        R_in = 3 * R_s
        R_out = 1000 * R_s

        # Integrate Planck emission over disk radius
        r_vals = np.logspace(np.log10(R_in), np.log10(R_out), 100)
        integrand = lambda r: 2 * np.pi * r * np.cos(self.theta) * \
                    self.planck_function(nu * (1+0.0), self.disk_temperature(r, M_bh, M_dot))
        # We integrate in log space for better accuracy:
        vals = [integrand(r) for r in r_vals]
        integral = np.trapz(np.array(vals) * r_vals, x=np.log(r_vals))

        # Divide by distance^2
        d_L = self.luminosity_distance(self.z)
        F_nu = integral / d_L**2

        # Apply extinction
        atten = self.extinction_factor(nu, ebv)
        return F_nu * atten

    def simulate_spectrum(self, freqs, M_bh, M_dot, ebv):
        """
        Generate model spectrum (nu * F_nu) at observer frequencies.
        freqs: array of frequencies [Hz] in the **observer** frame.
        """
        # Shift to emitted frame by (1+z), compute F_nu, then multiply by nu(1+z)
        # to compare directly with νF_ν data.
        spectrum = np.array([self.flux_density(nu*(1+self.z), M_bh, M_dot, ebv) for nu in freqs* (1 + self.z)])
        return spectrum

  def log_prior(self, params):
        """
        Prior on [log10(M_bh [g]), log10(Mdot [g/s]), log_f, E(B-V)].
        """
        logM, logMdot, logf, ebv = params
        M_linear = 10**logM
        Mdot_linear = 10**logMdot

        Mdot_Edd = (L_Edd_factor * (M_linear/M_sun)) / (0.1 * c**2)
        logM_min, logM_max = self.M_range
        logMdot_min = self.Mdot_min
        ebv_min, ebv_max = self.ebv_range

        if (logM_min < logM < logM_max and
            logMdot_min < logMdot < np.log10(Mdot_Edd) and
            self.logf_range[0] < logf < self.logf_range[1] and
            ebv_min <= ebv <= ebv_max):  # <- Prior now respects custom bounds
            return 0.0
        return -np.inf

    def log_likelihood(self, params):
        """
        Gaussian log-likelihood comparing log10(data) to log10(model).
        params = [logM, logMdot, logf, ebv].
        """
        logM, logMdot, logf, ebv = params
        M_linear = 10**logM
        Mdot_linear = 10**logMdot

        # Model fluxes at observed frequencies (ν F_ν)
        model_flux = self.simulate_spectrum(self.freq, M_linear, Mdot_linear, ebv)*(10**self.freq_log)* (1 + self.z)
        model_log = np.log10(model_flux)  # avoid log(0)

        # Combine errors: data error plus additional (f-factor) noise in log space
        sigma2 = self.flux_log_err**2 + np.exp(2 * logf)
        chi2 = np.sum((self.flux_log - model_log)**2 / sigma2 + np.log(sigma2))
        return -0.5 * chi2

    def log_posterior(self, params):
        """Sum of log-prior and log-likelihood."""
        lp = self.log_prior(params)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(params)

    def run_sampler(self):
        """Run the MCMC sampler (emcee) and store the flat chain."""
        ndim = 4
        # Initial positions of walkers: add small random scatter around initial guess
        init_log = [np.log10(self.initial_values[0]*M_sun),
                    np.log10(self.initial_values[1]*M_sun/(365*24*3600)),
                    self.initial_values[2],
                    self.initial_values[3]]
        pos = init_log + 1e-2 * np.random.randn(self.nwalkers, ndim)

        # Set up HDF5 backend
        out_dir = "./"
        filename = os.path.join(out_dir, "output_ext.h5")
        backend = emcee.backends.HDFBackend(filename)
        backend.reset(self.nwalkers, ndim)

        # Run sampler with multiprocessing
        with Pool() as pool:
            sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.log_posterior,
                                            args=(), pool=pool, backend=backend)
            sampler.run_mcmc(pos, self.nsteps, progress=True)

        # Flatten the chain (discard burn-in and thin)
        self.samples = sampler.get_chain(discard=self.nburn, thin=15, flat=True)

    def plot_corner(self):
        """Make a corner plot of the posterior samples."""
        import corner
        # Convert to physical units for plotting:
        samples = self.samples.copy()
        samples[:,0] = 10**samples[:,0] / M_sun
        samples[:,1] = 10**samples[:,1] * (365*24*3600) / M_sun
        labels = ["$M_{BH}$ ($M_\\odot$)", "$\\dot{M}$ ($M_\\odot$/yr)",
                  "log$_{10}(f)$", "E(B-V)"]
        fig = corner.corner(samples, labels=labels, truths=self.initial_values)
        fig.savefig("corner_plot_ext.png")
        plt.close(fig)

    def plot_overlay(self):
        """Overlay plot of random model spectra on top of data."""
        freqs = self.frequencies
        inds = np.random.choice(len(self.samples), size=self.overlay_number, replace=False)

        plt.figure(figsize=(8,6))
        for ind in inds:
            logM, logMdot, logf, ebv = self.samples[ind]
            M = 10**logM
            Mdot = 10**logMdot
            # Compute spectrum (nu*F_nu) and plot
            spec = self.simulate_spectrum(freqs, M, Mdot, ebv)
            plt.loglog(freqs, spec, color='orange', alpha=0.3)
        # Plot observed data points
        plt.errorbar(self.freq, self.flux, yerr=self.flux_err, fmt='o', color='blue',
                     label="Observed Data")
        plt.xlabel("Frequency [Hz]")
        plt.ylabel(r"$\nu F_\nu$ [erg s$^{-1}$ cm$^{-2}$]")
        plt.title("Model Spectra overlaid on Data (with Dust Extinction)")
        plt.legend()
        plt.grid(which="both", ls="--")
        plt.savefig("model_overlay_ext.png")
        plt.close()

    def summarize_posteriors(self):
        """
        Compute median and ±1σ for each parameter, output to CSV and print.
        Columns: M_bh, Mdot, log10(f), E(B-V).
        """
        # Copy samples and convert first two to physical units (M_sun, M_sun/yr)
        flat = self.samples.copy()
        flat[:,0] = 10**flat[:,0] / M_sun
        flat[:,1] = 10**flat[:,1] * (365*24*3600) / M_sun

        summary = {}
        names = ["M_bh", "Mdot", "logf", "ebv"]
        for i, name in enumerate(names):
            p16, p50, p84 = np.percentile(flat[:, i], [16, 50, 84])
            summary[name] = {"median": p50,
                             "-1sigma": p50 - p16,
                             "+1sigma": p84 - p50}
        df = pd.DataFrame(summary).T
        df.to_csv("fit_results_ext.csv")
        print(df)
        return df
