"""
cluster_de_density_model.py

Run: in terminal: python cluster_de_density_model.py
Or run cells in Jupyter / Colab.

Dependencies:
- numpy
- matplotlib
- scipy (for nicer color maps / stats; optional)
"""

import numpy as np
import matplotlib.pyplot as plt

# ------------------ Constants ------------------
G = 6.674e-11            # m^3 kg^-1 s^-2
kB = 1.381e-23           # J K^-1
c = 3e8                  # m s^-1
Msun = 1.989e30          # kg
Mpc = 3.086e22           # m
mp = 1.673e-27           # kg (proton mass)
mu_default = 0.6         # mean molecular weight (typical intracluster gas)
Lambda_obs = 1.11e-52    # m^-2 (approx observed cosmological constant)

# ------------------ Core function ------------------
def T_vir(M, R, Lambda, mu=mu_default):
    """
    Virial temperature (K) including a Λ-term correction.
    Inputs:
      M: mass in kg (can be scalar or numpy array)
      R: radius in m (same shape as M or scalar)
      Lambda: cosmological constant in m^-2 (same shape or scalar)
      mu: mean molecular weight (dimensionless)
    Returns:
      T in Kelvin (same shape as inputs)
    """
    term_grav = (G * M) / R
    term_DE = (1/6) * Lambda * (c**2) * (R**2)
    T = (mu * mp / (3 * kB)) * (term_grav - term_DE)
    return T

# ------------------ Scenario 1: Λ Sweep ------------------
def scenario_lambda_sweep_dynamic():
    # Range of cluster masses (in solar masses)
    M_values = [1e14, 3e14, 1e15]          # small, medium, large clusters
    R_values = [1.0, 2.0, 3.0]             # corresponding radii (Mpc)
    Lambda_values = np.linspace(0, 2.0e-52, 400)

    plt.figure(figsize=(8,6))

    for M_sun, R_Mpc in zip(M_values, R_values):
        M = M_sun * Msun
        R = R_Mpc * Mpc
        T_vals = T_vir(M, R, Lambda_values)
        plt.plot(Lambda_values * 1e52, T_vals/1e7, lw=2,
                 label=f"M={M_sun:.0e} M☉, R={R_Mpc:.1f} Mpc")

    plt.axvline(Lambda_obs*1e52, color='black', ls='--', label='Λ_obs')
    plt.xlabel("Λ (×10⁻⁵² m⁻²)")
    plt.ylabel("T_vir (×10⁷ K)")
    plt.title("Dynamic Λ Sweep — Multiple Cluster Sizes")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

# ------------------ Scenario 2: Mass–Radius Grid (Heatmap) ------------------
def scenario_mass_radius_grid():
    M_vals = np.array([1e14, 3e14, 5e14, 1e15]) * Msun
    R_vals = np.array([0.8, 1.5, 2.5, 4.0]) * Mpc
    frac_pct = np.zeros((len(M_vals), len(R_vals)))

    for i, M in enumerate(M_vals):
        for j, R in enumerate(R_vals):
            T0 = T_vir(M, R, 0)
            Tlam = T_vir(M, R, Lambda_obs)
            frac_pct[i, j] = (T0 - Tlam) / T0 * 100.0

    # plot heatmap
    plt.figure(figsize=(7,5))
    im = plt.imshow(frac_pct, origin='lower', cmap='plasma',
                    extent=[R_vals[0]/Mpc, R_vals[-1]/Mpc, M_vals[0]/Msun/1e14, M_vals[-1]/Msun/1e14],
                    aspect='auto')
    cbar = plt.colorbar(im)
    cbar.set_label('% Decrease in T due to DE')
    plt.xlabel('Radius (Mpc)')
    plt.ylabel('Mass (10^14 M☉)')
    plt.title('Scenario 2 — Fractional DE Effect on T (Λ = Λ_obs)')
    plt.tight_layout()
    plt.show()

# ------------------ Scenario 3: Mass Scaling with R ∝ M^(1/3) ------------------
def scenario_mass_scaling():
    M_range = np.logspace(14, 16, 60) * Msun  # 1e14 -> 1e16 M☉
    R0 = 2.0 * Mpc
    M_ref = 5e14 * Msun
    # Simple scaling: R(M) = R0 * (M/M_ref)^(1/3)
    R_range = R0 * (M_range / M_ref)**(1.0/3.0)
    T_noDE = T_vir(M_range, R_range, 0)
    T_withDE = T_vir(M_range, R_range, Lambda_obs)

    plt.figure(figsize=(7,5))
    plt.loglog(M_range / Msun, T_noDE / 1e7, label='Λ = 0')
    plt.loglog(M_range / Msun, T_withDE / 1e7, label='Λ = Λ_obs')
    plt.xlabel('Cluster Mass (M☉)')
    plt.ylabel('T_vir (×10⁷ K)')
    plt.title('Scenario 3 — Mass Scaling (R ∝ M^{1/3})')
    plt.legend()
    plt.grid(True, which='both', alpha=0.3)
    plt.tight_layout()
    plt.show()

# ------------------ Scenario 4: Redshift Evolution (R ∝ 1/(1+z)) ------------------
def scenario_redshift_evolution():
    z = np.linspace(0.0, 3.0, 300)
    M = 5e14 * Msun
    R0 = 2.0 * Mpc

    def R_of_z(R0, z):
        return R0 / (1.0 + z)   # simplified scaling for illustration

    R_z = R_of_z(R0, z)
    T_noDE = T_vir(M, R_z, 0)
    T_withDE = T_vir(M, R_z, Lambda_obs)

    plt.figure(figsize=(7,5))
    plt.plot(z, T_noDE / 1e7, label='Λ = 0')
    plt.plot(z, T_withDE / 1e7, label='Λ = Λ_obs')
    plt.xlabel('Redshift z')
    plt.ylabel('T_vir (×10⁷ K)')
    plt.title('Scenario 4 — Virial Temperature Evolution with Redshift')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

# ------------------ Scenario 5: Monte Carlo Uncertainty Propagation ------------------
def scenario_monte_carlo(n_samples=10000):
    np.random.seed(42)
    M0 = 5e14 * Msun
    R0 = 2.0 * Mpc
    sigma_frac_M = 0.2   # 20% uncertainty
    sigma_frac_R = 0.2

    M_samples = np.random.normal(M0, sigma_frac_M * M0, size=n_samples)
    R_samples = np.random.normal(R0, sigma_frac_R * R0, size=n_samples)

    T_noDE = T_vir(M_samples, R_samples, 0)
    T_withDE = T_vir(M_samples, R_samples, Lambda_obs)

    # quick stats
    mean0, std0 = np.mean(T_noDE), np.std(T_noDE)
    mean1, std1 = np.mean(T_withDE), np.std(T_withDE)
    print("Monte Carlo (n={} samples)".format(n_samples))
    print(f"Λ=0: mean T = {mean0:.3e} K, std = {std0:.3e} K")
    print(f"Λ=Λ_obs: mean T = {mean1:.3e} K, std = {std1:.3e} K")
    print("Mean fractional shift due to DE: {:.3e}%".format((mean0 - mean1)/mean0 * 100.0))

    # plot histograms
    plt.figure(figsize=(8,5))
    plt.hist(T_noDE/1e7, bins=60, alpha=0.6, label='Λ=0')
    plt.hist(T_withDE/1e7, bins=60, alpha=0.6, label='Λ=Λ_obs')
    plt.xlabel('T_vir (×10⁷ K)')
    plt.ylabel('Count')
    plt.title('Scenario 5 — Monte Carlo: Distribution of T_vir with Measurement Uncertainties')
    plt.legend()
    plt.tight_layout()
    plt.show()

# ------------------ Optional: Scenario 6: varying w (dynamic DE) ------------------
def scenario_w_variation():
    z = np.linspace(0.0, 3.0, 300)
    M = 5e14 * Msun
    R0 = 2.0 * Mpc

    def R_of_z(R0, z):
        return R0 / (1.0 + z)

    def Lambda_w(Lambda0, z, w):
        # Simplified scaling for DE density: rho_DE(z) = rho_DE0 * (1+z)^{3(1+w)}
        # Translate to effective Lambda_eff(z) ~ Lambda0 * (1+z)^{3(1+w)}
        return Lambda0 * (1.0 + z)**(3.0 * (1.0 + w))

    for w in [-1.1, -1.0, -0.9]:
        Lambda_eff = Lambda_w(Lambda_obs, z, w)
        T_vals = T_vir(M, R_of_z(R0, z), Lambda_eff)
        plt.plot(z, T_vals / 1e7, label=f'w={w}')

    plt.xlabel('Redshift z')
    plt.ylabel('T_vir (×10⁷ K)')
    plt.title('Scenario 6 — Effect of DE equation-of-state (w) on T_vir(z)')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

# ------------------ Main runner ------------------
def main():
    print("Running cluster DE-virial simulation suite...")
    scenario_lambda_sweep_dynamic()
    scenario_mass_radius_grid()
    scenario_mass_scaling()
    scenario_redshift_evolution()
    scenario_monte_carlo(n_samples=5000)    # reduce samples if slow
    scenario_w_variation()
    print("Done.")

if __name__ == "__main__":
    main()
