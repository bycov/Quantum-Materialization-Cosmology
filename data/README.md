# Data Sources

This folder contains the datasets used in the QMG analysis.

## `Pantheon+SH0ES.dat`
- **Source:** Pantheon+ Supernovae compilation
- **Reference:** Brout et al. (2022), ApJ 938, 110 [arXiv:2202.04077]
- **Content:** Contains redshifts, distance moduli and errors for Type Ia supernovae.

## External Data Used in Code (Not Stored Here)

The following data are not stored as files but are hardcoded in `../code/QMG_MCMC.py` as published numerical results:

- **DESI BAO:** 7 measurements of `D_V(z)/r_d` from Adame et al. (2024) [arXiv:2404.03002] [citation:8]
- **fσ₈(z):** 15 measurements from 6dFGS, SDSS, BOSS, eBOSS, WiggleZ compiled by Alam et al. (2017) [arXiv:1607.03155]
- **KiDS-1000 (`S_8`):** `S_8 = 0.766 ± 0.020` from Heymans et al. (2021) [arXiv:2007.15632] [citation:2][citation:7]
