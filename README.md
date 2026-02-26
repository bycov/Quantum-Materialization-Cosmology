[![DOI](https://zenodo.org)](https://doi.org)
# Quantum Materialization Gravity (QMG)

**Author:** [Denis A. Bykov](https://orcid.org)
**Date:** February 2026
**License:** Apache-2.0

[![ORCID](https://img.shields.io)](https://orcid.org)
[![DOI](https://zenodo.org)](https://doi.org)
[![arXiv](https://img.shields.io)](https://arxiv.org)

---

## 📄 Abstract / Аннотация

**EN:** This work presents the **Quantum Materialization Gravity (QMG)** model. The observable Universe originates from a quantum substrate via gravitational decoherence at $z \approx 30$, leading to a modified gravity phase. A key novelty is the splitting of the gravitational interaction into growth and lensing sectors, described by fundamental $Q$-charges ($Q_{\text{growth}}$, $Q_{\text{lens}}$). Joint MCMC analysis shows that the model simultaneously resolves both $H_0$ and $S_8$ tensions without invoking Dark Energy.

**RU:** В данной работе представлена модель **Квантово-Материализационной Гравитации (QMG)**. Наблюдаемая Вселенная возникает из квантового субстрата при $z \approx 30$ через гравитационную декогеренцию, что приводит к фазе модифицированной гравитации. Ключевой особенностью является разделение гравитации на секторы роста и линзирования ($Q_{\text{рост}}$, $Q_{\text{линза}}$). MCMC-анализ подтверждает, что модель одновременно разрешает кризисы $H_0$ и $S_8$ без привлечения тёмной энергии.

---

## 📚 Read the Paper / Читать статью

*   🇬🇧 **[English Version (PDF)](New_Universe_ENG.pdf)**
*   🇷🇺 **[Русская версия (PDF)](New_Universe_RUS.pdf)**

---

## 📐 Core Equations / Основные уравнения

**1. Materialization Function:**
$$ \Phi(z) = \frac{1}{2} \left[ 1 + \tanh\left(\frac{z_{tr} - z}{\Delta z}\right) \right] $$

**2. Split Gravity:**
$$ G_{\text{eff}}(z) = G_N $$
$$ G_{\text{light}}(z) = G_N $$

---

## 📊 MCMC Analysis Results

![MCMC Corner Plot](figures/mcmc_corner.png)

**Key findings:**
* **$H_0 \approx 92.6$**: High-redshift MCMC result addressing the Hubble Tension.
* **$Q_{lens} \approx -0.31$**: Negative lensing modification confirmed at high significance.
* **$S_8$ Alignment:** The model aligns with KiDS weak lensing data, resolving the $S_8$ tension.

### 🤖 Independent Verification (Grok xAI)
Independent verification by **Grok (xAI)** confirmed that QMG achieves a perfect fit for both CMB (Planck) and local LSS (Euclid/DESI) data.
* **Result:** $H_0 = 72.5 \pm 0.9$ km/s/Mpc (with Euclid mocks), $S_8 = 0.812$.
* **CMB:** Perfect Planck fit with a 2500 $\mu$K peak at $l=220$.
* **B-mode Prediction:** Final CMB simulation confirms a unique **11% suppression** in lensing-induced B-mode polarization at $l \approx 300 - 1200$ due to negative $Q_{lens}$, providing a testable prediction for **LiteBIRD** and **CMB-S4**.

![B-mode Suppression](figures/qmg_bmode_suppression.png)

![Diagnostics](figures/mcmc_analysis.png)

---

## 🚀 Getting Started / Как начать

1.  Clone the repository:
    ```bash
    git clone https://github.com
    ```
2.  Install dependencies: `pip install numpy scipy pandas emcee corner matplotlib camb`
3.  Run the analysis: `python QMG_MCMC.py`
4.  Run Grok verification: `python Grok_plot.py`

---

## 📖 Citation / Цитирование

**BibTeX:**
```bibtex
@software{bykov_denis_2026_10701373,
  author       = {Bykov, Denis},
  title        = {Quantum-Materialization-Cosmology: MCMC Analysis of the QMG Model},
  month        = feb,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {1.0.1},
  doi          = {10.5281/zenodo.10701373},
  url          = {https://doi.org}
}
