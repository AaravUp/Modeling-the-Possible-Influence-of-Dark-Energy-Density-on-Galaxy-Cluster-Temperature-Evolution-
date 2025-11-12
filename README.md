# Modeling-the-Possible-Influence-of-Dark-Energy-Density-on-Galaxy-Cluster-Temperature-Evolution-
Simulation and analysis of dark energy effects on galaxy cluster thermodynamics.
⚙️ How to Use
1. Install Dependencies

Make sure you have Python ≥ 3.9 installed.

Install the required libraries:

pip install numpy matplotlib pandas

2. Run the Simulations

To execute all scenarios at once:

python cluster_de_density_model.py


Each scenario will automatically generate a plot window.
You can save figures by adding plt.savefig("filename.png", dpi=300) inside any scenario function.

🌌 Simulation Scenarios
Scenario 1 — Λ Sweep

Goal: Observe how virial temperature changes with increasing dark energy density (Λ).
Description:

Keeps cluster mass and radius fixed.

Varies Λ from 0 → 2×10⁻⁵² m⁻².

Demonstrates that Λ introduces a small “cooling” correction.
Output: Temperature vs Λ plot (Figure 1).

Scenario 2 — Mass–Radius Grid

Goal: Quantify the percentage decrease in temperature across different cluster sizes.
Description:

Creates a 2D grid of mass (10¹⁴–10¹⁵ M☉) and radius (1–4 Mpc).

Calculates fractional ΔT between Λ=0 and Λ=Λₒᵦₛ.
Output: Heatmap showing which clusters are most sensitive to dark energy (Figure 2).

Scenario 3 — Mass Scaling

Goal: Examine the scaling T ∝ M ^ (2/3) with and without Λ.

Description:

Relates radius to mass using R∝ M ^ (1/3).

Plots log–log relation of temperature vs mass.

Verifies that Λ slightly lowers normalization but doesn’t change slope.

Output: Log–log T–M plot (Figure 3).

Scenario 4 — Redshift Evolution

Goal: Model how cluster temperature evolves with cosmic time.

Description:

Contracts cluster size with redshift (R∝1/(1+z)).

Computes T(z) for Λ=0 and Λ=Λₒᵦₛ.

Shows that dark energy flattens the T–z curve at low redshift.
Output: Temperature vs Redshift plot (Figure 4).

Scenario 5 — Monte Carlo Uncertainty Propagation

Goal: Test if the Λ effect is observable given realistic measurement noise.

Description:

Introduces ±20% random errors in M and R.

Runs 10,000 samples for Λ=0 and Λ=Λₒᵦₛ.

Compares statistical distributions of T.

Output: Histograms showing that Λ effects are orders of magnitude smaller than noise (Figure 5).

Scenario 6 — Observational Comparison

Goal: Compare modeled and observed cluster temperatures.
Description:

Loads real data (from vizier_votable.tsv or cluster_data.csv).

Computes modeled virial temperature using T_vir() function.

Plots observed vs modeled temperatures and correlation.

Output: Scatter plot comparing data and simulation (Figure 6).
