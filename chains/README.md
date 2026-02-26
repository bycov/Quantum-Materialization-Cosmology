# MCMC Production Chains (QMG Model)

This directory contains the final MCMC output for the Quantum Materialization Gravity (QMG) analysis.

### Files Description:
* **chains_QMG.npy** — Raw MCMC chains (NumPy binary format). Includes 300 steps for 32 walkers across 6 parameters ($H_0, \Omega_m, z_{tr}, Q_{growth}, Q_{lens}, \alpha$).
* **chains_QMG_labels.txt** — Plain text file mapping columns in the `.npy` file to physical parameters.

### Chain Statistics:
* **H0** ≈ 92.6 km/s/Mpc
* **Q_lens** ≈ -0.31 (Negative lensing confirmed at >2σ)
* **Status:** Converged (Burn-in removed)

These chains are ready for cross-correlation with CMB simulations.
