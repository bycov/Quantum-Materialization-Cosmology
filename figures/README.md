# MCMC Analysis Figures (QMG Model)

This folder contains the visualization of results from the Markov Chain Monte Carlo (MCMC) analysis for the Quantum Materialization Gravity model.

### Main Diagnostic Plots
* **mcmc_corner.png** — Corner plot showing posterior distributions and correlations between $H_0$, $\Omega_m$, $z_{tr}$, $Q_{growth}$, $Q_{lens}$, and $\alpha$.
* **mcmc_analysis.png** — Dashboard including $S_8$ distribution, $\Omega_m$ vs $Q_{lens}$ correlation, $f\sigma_8$ growth data fit, and MCMC trace plots.

### Individual Plots (for Paper/Presentations)
* **fs8_fit.png** — Growth of structure data ($f\sigma_8$) vs QMG model prediction.
* **bao_fit.png** — Baryon Acoustic Oscillations data vs model distance predictions.
* **s8_distribution.png** — Posterior distribution of the $S_8$ parameter compared to Planck and KiDS-1000.
* **corner.png** — Legacy version of the corner plot.

---
**Note:** All plots were generated using `QMG_MCMC.py` with 300 steps and 32 walkers.
