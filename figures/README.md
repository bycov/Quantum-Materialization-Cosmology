# MCMC Analysis Figures (QMG Model)

This folder contains visualization of results from the Markov Chain Monte Carlo (MCMC) analysis for the Quantum Materialization Gravity model.

## Main Diagnostic Plots

- **`mcmc_corner.png`** — Corner plot showing posterior distributions and correlations between parameters: $H_0$, $\Omega_m$, $z_{tr}$, $Q_{\text{growth}}$, $Q_{\text{lens}}$, and $\alpha$.
- **`mcmc_analysis.png`** — Dashboard including $S_8$ distribution, $\Omega_m$ vs $Q_{\text{lens}}$ correlation, $f\sigma_8(z)$ growth fit, and MCMC trace plots.
- **`qmg_bmode_suppression.png`** — **Key prediction:** 11% suppression in lensing B-modes at $l \approx 300-1200$ due to negative $Q_{\text{lens}} = -0.16$.

## Data Comparison Plots

- **`fs8_fit.png`** — Structure growth: QMG model (blue) vs $f\sigma_8(z)$ data (red points). Perfect plateau at $f\sigma_8 \approx 0.45$.
- **`bao_fit.png`** — BAO distance: DESI data vs QMG prediction.
- **`s8_distribution.png`** — $S_8$ parameter distribution: QMG (blue histogram) compared to Planck (red) and KiDS-1000 (green) constraints.

## Notes

- All plots generated with `code/QMG_MCMC.py`
- MCMC chains available in `../chains/`
- Mock predictions for future experiments in `../predictions/`
