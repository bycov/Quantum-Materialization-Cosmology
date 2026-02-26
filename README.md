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

**EN:** This work presents the **Quantum Materialization Gravity (QMG)** model, an extension of the initial QMC framework. The observable Universe originates from a quantum substrate via gravitational decoherence at $z \approx 30$, leading to a modified gravity phase. A key novelty is the splitting of the gravitational interaction into growth and lensing sectors, described by fundamental $Q$-charges ($Q_{\text{growth}}$, $Q_{\text{lens}}$). A joint MCMC analysis of **DESI BAO**, **Pantheon+**, $f\sigma_8(z)$, and **KiDS-1000** data yields excellent agreement. The model simultaneously resolves both the $H_0$ and $S_8$ tensions without invoking Dark Energy. A crucial finding is the negative lensing charge, $Q_{\text{lens}} \approx -0.32$, indicating a weaker gravitational effect on light during the materialization epoch.

**RU:** В данной работе представлена модель **Квантово-Материализационной Гравитации (QMG)** — расширение первоначальной модели QMC. Наблюдаемая Вселенная возникает из квантового субстрата при $z \approx 30$ через гравитационную декогеренцию, что приводит к фазе модифицированной гравитации. Ключевой особенностью является разделение гравитационного взаимодействия на секторы роста и линзирования, описываемые фундаментальными $Q$-зарядами ($Q_{\text{рост}}$, $Q_{\text{линза}}$). Совместный MCMC-анализ данных **DESI BAO**, **Pantheon+**, $f\sigma_8(z)$ и **KiDS-1000** демонстрирует отличное согласие. Модель одновременно разрешает оба космологических кризиса ($H_0$ и $S_8$ tensions) без привлечения тёмной энергии. Ключевой результат — отрицательный линзирующий заряд $Q_{\text{линза}} \approx -0.32$, указывающий на ослабление гравитации для света в эпоху материализации.

---

## 📚 Read the Paper / Читать статью

*   🇬🇧 **[English Version (PDF)](New_Universe_ENG.pdf)**
*   🇷🇺 **[Русская версия (PDF)](New_Universe_RUS.pdf)**

---

## 📐 Core Equations / Основные уравнения

**1. Materialization Function / Функция материализации:**
$$ \Phi(z) = \frac{1}{2} \left[ 1 + \tanh\left(\frac{z_{tr} - z}{\Delta z}\right) \right] $$

**2. Split Gravity / Разделение гравитации:**
$$ G_{\text{eff}}(z) = G_N $$
$$ G_{\text{light}}(z) = G_N $$

---

## 📊 MCMC Analysis Results

The QMG model parameters were constrained using Pantheon+ SNe Ia, BAO, and $f\sigma_8$ data.

![MCMC Corner Plot](figures/mcmc_corner.png)

**Key findings:**
* **$H_0 \approx 96.5$**: Significantly higher than the standard $\Lambda$CDM value, addressing the Hubble Tension.
* **$Q_{lens} \approx -0.32$**: Strongly suggests a negative lensing modification (3$\sigma$ significance).
* **$S_8$ Alignment:** The model aligns closely with KiDS weak lensing data, showing a preference for a lower $S_8$ compared to the Planck CMB baseline.

### Model Diagnostics & Tensions
![Diagnostics](figures/mcmc_analysis.png)

* **$H_0$ Convergence:** Trace plots indicate stable MCMC convergence after ~150 steps.
* **Structure Growth:** Significant deviation in $f\sigma_8$ at low redshifts, providing a unique signature of the QMG model.
* **Statistics:** $\chi^2_{\text{red}}(f\sigma_8) = 1.01$

---

## 📁 Repository Content / Состав репозитория


| File / Файл | Description / Описание |
| :--- | :--- |
| `QMG_MCMC.py` | Main MCMC analysis code / Основной код анализа |
| `sigma8_normalization.py` | $\sigma_8$ normalization via CAMB / Нормировка через CAMB |
| `Pantheon+SH0ES.dat` | Supernovae dataset / Данные сверхновых |
| `figures/` | Folder containing all plots / Папка с графиками |

---

## 🚀 Getting Started / Как начать

1.  Clone the repository:
    ```bash
    git clone https://github.com
    ```
2.  Install dependencies: `pip install numpy scipy pandas emcee corner matplotlib`
3.  Run the analysis: `python QMG_MCMC.py`

---

## 📖 Citation / Цитирование

If you use this code or model in your research, please cite it as follows:

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
