# Quantum-Materialization-Cosmology (QMC)

**Author:** D.A. Bykov  
**Date:** February 2026  
**License:** Apache-2.0

---

## 📄 Abstract / Аннотация

**EN:** This paper proposes the **Quantum Materialization Cosmology (QMC)** model, where the observable Universe originates from a quantum substrate through gravitational-induced decoherence (materialization) at $z \approx 40$, rather than a Big Bang singularity. In this framework, the materialization process leads to a **Modified Gravity** phase, effectively explaining the late-time accelerated expansion of the Universe without invoking Dark Energy. MCMC analysis using **DESI BAO**, **Pantheon+**, and $f\sigma_8(z)$ datasets demonstrates excellent statistical agreement with observations ($\chi^2/dof = 1.01$). The model yields a Hubble constant of $H_0 = 82.9 \pm 9.7$ km/s/Mpc, providing a potential pathway to mitigate the $H_0$ tension.

**RU:** В данной работе предлагается модель **Квантовой материализации Вселенной (QMC)**, в которой наблюдаемая Вселенная возникает из квантового субстрата в результате гравитационно-индуцированной декогеренции (материализации) при $z \approx 40$, заменяя сингулярность Большого взрыва. Процесс материализации приводит к фазе **Модифицированной гравитации**, эффективно объясняя ускоренное расширение Вселенной без привлечения Тёмной энергии. MCMC-анализ данных **DESI BAO**, **Pantheon+** и $f\sigma_8(z)$ демонстрирует отличное статистическое согласие с наблюдениями ($\chi^2/dof = 1.01$). Модель дает значение постоянной Хаббла $H_0 = 82.9 \pm 9.7$ км/с/Мпк, что открывает путь к решению проблемы $H_0$ tension.

---

## 📚 Read the Paper / Читать статью

*   🇬🇧 **[English Version (PDF)](New_Universe_ENG.pdf)** — Full theoretical paper in English.
*   🇷🇺 **[Русская версия (PDF)](New_Universe_RUS.pdf)** — Полный текст научной работы на русском языке.

---

## 📐 Key Equations / Основные уравнения

The evolution of the effective gravitational constant $G_{eff}(z)$ is governed by the materialization function $\Phi(z)$:

1. **Transition Function / Функция перехода:**
$$\Phi(z) = \frac{1}{2} \left[ 1 + \tanh\left(\frac{z_{tr} - z}{\Delta z}\right) \right]$$

2. **Modified Gravity / Модифицированная гравитация:**
$$G_{eff}(z) = G_N [1 + \beta \cdot \Phi(z)]$$

---

## 📊 MCMC Analysis Results / Результаты анализа


| Parameter | Value (68% CL) |
| :--- | :--- |
| **$H_0$** | $82.9 \pm 9.7$ km/s/Mpc |
| **$\Omega_m$** | $0.30 \pm 0.14$ |
| **$z_{tr}$** | $39.8 \pm 13.4$ |
| **$\beta$** | $0.76 \pm 0.50$ |

**Statistics:** $\chi^2/dof = 1.01$

---

## 📁 Repository Content / Состав репозитория

*   `New_Universe_ENG.pdf` — Paper (English).
*   `New_Universe_RUS.pdf` — Статья (Русский).
*   `Modified_gravity.py` — Python script for MCMC calculations.
*   `Pantheon+SH0ES.dat` — Observational dataset.

---
*Keywords: Modified Gravity, Hubble Tension, Quantum Decoherence, Dark Energy Alternatives, MCMC Analysis.*
