# Nested Mixed Model Analysis in R — Rats Glycogen Experiment

> **Author:** Ms. Mrunalini (Data Science Trainer)  
> **Contact:** mrunalini0107@gmail.com | Mumbai – 400095

---

## 📋 Table of Contents

1. [Project Description](#project-description)
2. [Experiment Overview](#experiment-overview)
3. [Learning Objectives](#learning-objectives)
4. [Project Structure](#project-structure)
5. [Requirements](#requirements)
6. [Key Concepts Covered](#key-concepts-covered)
7. [Workflow Summary](#workflow-summary)
8. [Results & Conclusion](#results--conclusion)
9. [How to Run](#how-to-run)

---

## Project Description

This project demonstrates a **complete statistical workflow** for analysing hierarchically structured biological data using **Linear Mixed-Effects Models (LME)** in R. The dataset comes from a classic nested experimental design in which glycogen concentration is measured across multiple rats, liver parts, and treatment groups.

The notebook walks through every stage of a rigorous mixed-model analysis: data inspection, correct model specification, assumption diagnostics, random and fixed effects interpretation, model selection using AIC and likelihood ratio tests, post-hoc comparisons with estimated marginal means, and simulation-based power analysis.

A key pedagogical feature of the project is the **side-by-side comparison of correct and incorrect model specifications**, making it an excellent reference for anyone learning nested random-effects structures in `lme4`.

---

## Experiment Overview

| Element | Detail |
|---|---|
| **Response variable** | Glycogen concentration |
| **Fixed effect** | Treatment (3 levels: Treat1, Treat2, Treat3) |
| **Random effect 1** | Rat nested within Treatment |
| **Random effect 2** | Liver part (Left, Middle, Right) nested within Rat |
| **Replication** | 2 measurements per liver part |
| **Total observations** | 3 treatments × 2 rats × 3 liver parts × 2 replicates = **36** |

### Hierarchical Structure

```
Treatment (Fixed)
 └── Rat (Random, nested in Treatment)
       └── Liver Part (Random, nested in Rat)
             └── Measurement replicate (Residual)
```

### Mathematical Model

$$Y_{ijkl} = \mu + \alpha_i + b_{j(i)} + c_{k(j(i))} + \varepsilon_{ijkl}$$

Where:
- $\mu$ = overall mean glycogen level  
- $\alpha_i$ = fixed effect of Treatment $i$  
- $b_{j(i)} \sim N(0, \sigma^2_{Rat})$ = random effect of Rat nested within Treatment  
- $c_{k(j(i))} \sim N(0, \sigma^2_{Liver})$ = random effect of Liver nested within Rat  
- $\varepsilon_{ijkl} \sim N(0, \sigma^2)$ = residual error  

---

## Learning Objectives

By working through this notebook, you will be able to:

1. **Understand nested vs. crossed experimental designs** and correctly identify the structure in real biological data.
2. **Specify linear mixed-effects models** in R using `lme4`, with proper fixed and random effect syntax for nested hierarchies.
3. **Diagnose model assumptions** including convergence, singularity, homoscedasticity (Levene's test), and normality of residuals (Shapiro–Wilk, Q–Q plots).
4. **Interpret variance components** from random effects and understand what rat-level vs. liver-level vs. residual variability means biologically.
5. **Perform model selection** using AIC and likelihood ratio tests (ANOVA), and understand when to use REML vs. ML estimation.
6. **Conduct and interpret post-hoc pairwise comparisons** using estimated marginal means (`emmeans`) with multiple comparison corrections (FDR, Tukey).
7. **Run simulation-based power analyses** for both fixed effects (Treatment) and random effects (Rat), and correctly interpret underpowered studies.
8. **Recognise and avoid common model mis-specification errors**, including treating a fixed factor as random and ignoring the correct nesting structure.
9. **Use diagnostic tools** from `DHARMa`, `performance`, `influence.ME`, and `car` to validate model adequacy.
10. **Communicate statistical findings clearly**, including when a non-significant result is inconclusive rather than evidence of no effect.

---

## Project Structure

```
R12_Nested_Mixed_Model_Analysis_-_Rats_Glycogen_Experiment.ipynb
│
├── 1. Nested Design Theory
│     ├── Nested vs. crossed designs
│     ├── Experiment description and factor table
│     └── Mathematical model formulation
│
├── 2. Data Inspection
│     ├── Raw data import (Rats2.xlsx)
│     ├── Factor conversion and unique rat ID creation
│     ├── Summary statistics
│     └── Visualisations (boxplots, nested structure plots)
│
├── 3. Model Specification
│     ├── Correct nested model: Treatment + (1|Treatment:Rat) + (1|Treatment:Rat:Liver)
│     ├── Wrong model 1: Treatment as random via (1|Treatment/Rat/Liver)
│     └── Wrong model 2: Non-unique Rat IDs with (1|Rat/Liver)
│
├── 4. Model Diagnostics
│     ├── Convergence and singularity checks (performance)
│     ├── Levene's test for homoscedasticity
│     ├── Q–Q plot and Shapiro–Wilk normality test
│     ├── Comprehensive check_model() diagnostics
│     ├── Cook's distance influence analysis (influence.ME)
│     └── DHARMa simulated residual diagnostics
│
├── 5. Summary Statistics & Inference
│     ├── lme4 model summary and variance components
│     ├── ANOVA table for fixed effects
│     ├── Estimated marginal means (emmeans)
│     ├── Pairwise comparisons (unadjusted, FDR, Tukey)
│     └── Random effects dotplots with confidence intervals
│
├── 6. Model Selection
│     ├── AIC comparison (REML and ML models)
│     ├── Likelihood ratio tests for random and fixed effects
│     └── Best model selection rationale
│
└── 7. Power Analysis
      ├── Simulation-based power for Rat random effect (Power = 0.34)
      └── Simulation-based power for Treatment fixed effect (Power = 0.68)
```

---

## Requirements

### R Packages

```r
install.packages(c(
  "readxl",       # Read Excel data
  "ggplot2",      # Data visualisation
  "lme4",         # Linear mixed-effects models
  "lmerTest",     # p-values for fixed effects
  "kableExtra",   # Formatted tables
  "performance",  # Model diagnostics
  "see",          # Visualisation companion to performance
  "car",          # Levene's test, Q–Q plots
  "emmeans",      # Estimated marginal means
  "DHARMa",       # Simulated residual diagnostics
  "influence.ME", # Cook's distance for mixed models
  "lattice"       # Random effects dotplots
))
```

### Data File

- `Rats2.xlsx` — Excel file containing glycogen measurements with columns: `Treatment`, `Rat`, `Liver`, `Glycogen`

---

## Key Concepts Covered

### Fixed vs. Random Effects
- **Fixed effects** (Treatment): Factors whose specific levels are of direct scientific interest. We want to compare Treat1, Treat2, and Treat3.
- **Random effects** (Rat, Liver): Factors whose levels are a random sample from a larger population. We want to account for, not compare, their variability.

### REML vs. Maximum Likelihood
| Setting | When to Use |
|---|---|
| `REML = TRUE` | Final model estimation; unbiased variance components |
| `REML = FALSE` | Comparing models with different **fixed** effects via AIC/LRT |

### Correct Nesting Syntax in lme4
```r
# Correct: explicit nesting with non-unique Rat IDs
lmer(Glycogen ~ Treatment + (1|Treatment:Rat) + (1|Treatment:Rat:Liver), data = Rats)

# Also correct: only when Rat IDs are globally unique
lmer(Glycogen ~ Treatment + (1|Rat) + (1|Rat:Liver), data = Rats)

# Wrong: Treatment treated as random effect
lmer(Glycogen ~ 1 + (1|Treatment/Rat/Liver), data = Rats)
```

---

## Workflow Summary

### Step 1 — Data Preparation
Import `Rats2.xlsx`, convert `Treatment`, `Rat`, and `Liver` to factors, and create a unique rat identifier (`Rats_Uniq`) to correctly represent six distinct rats rather than the repeated 1–2 coding.

### Step 2 — Exploratory Visualisation
Plot glycogen levels faceted by Treatment, Rat, and Liver to visually confirm the hierarchical nested structure before modelling.

### Step 3 — Model Fitting
Fit the primary nested mixed model with Treatment as a fixed effect and Rat-within-Treatment and Liver-within-Rat-within-Treatment as random intercepts.

### Step 4 — Diagnostics
Verify all assumptions: convergence (gradient check), no singularity, homoscedasticity (Levene, residuals vs fitted), normality (Shapiro–Wilk, Q–Q), absence of influential observations (Cook's distance < 0.5), and DHARMa residual uniformity.

### Step 5 — Inference
Examine variance components to understand sources of variability. Evaluate fixed effects via ANOVA and emmeans pairwise contrasts. Visualise random effect estimates with confidence intervals to assess their importance.

### Step 6 — Model Selection
Use AIC and likelihood ratio tests to compare full vs. reduced models. Confirm that including both Rat and Liver random effects (ModelR1) yields the lowest AIC, and that Treatment significantly improves model fit over a null model (p = 0.039).

### Step 7 — Power Analysis
Simulate 100 datasets to estimate power for detecting the Rat random effect (power = 0.34) and the Treatment fixed effect (power = 0.68). Both are below the recommended threshold of 0.80, indicating the study is underpowered.

### Step 8 — Model Comparison
Demonstrate and explain the consequences of three model specifications — correct, wrong (Treatment as random), and wrong (ignoring Treatment in nesting) — comparing REML values, variance estimates, and inferential validity.

---

## Results & Conclusion

### Variance Components (Best Model — ModelR1)

| Source | Variance | Std. Dev. | Interpretation |
|---|---|---|---|
| Rat within Treatment | 36.06 | 6.005 | Largest source — substantial biological variability between rats |
| Liver within Rat | 14.17 | 3.764 | Moderate variability across liver regions |
| Residual | 21.17 | 4.601 | Measurement and unexplained error |

### Fixed Effects (Treatment on Glycogen)

| Comparison | Estimate | p-value | Result |
|---|---|---|---|
| Treat2 vs Treat1 | +10.5 units | 0.213 | Not significant |
| Treat3 vs Treat1 | −5.3 units | 0.482 | Not significant |
| ANOVA (Treatment) | F = 2.93 | 0.197 | Not significant |

### Model Selection

The full model including both Rat and Liver random effects (ModelR1, AIC = 231.62) provided the best fit. Likelihood ratio testing confirmed that Treatment significantly improves the model over a null model (χ² = 6.50, df = 2, p = 0.039), while the addition of the Liver random effect beyond Rat alone was not statistically significant (p = 0.102), though it is retained on experimental design grounds.

### Power Analysis Summary

| Effect | Estimated Power | Interpretation |
|---|---|---|
| Rat random effect | 0.34 | Very low — high risk of missing true Rat variability |
| Treatment fixed effect | 0.68 | Low — study likely underpowered to detect true Treatment differences |

### Overall Conclusion

**Treatment does not have a statistically significant effect on glycogen concentration in this experiment** (F(2,3) = 2.93, p = 0.197). Observed mean differences between treatment groups — Treat2 slightly higher, Treat3 slightly lower than Treat1 — are not large enough relative to the substantial between-rat variability (σ²_Rat = 36.06) to reach significance.

However, this conclusion must be interpreted cautiously. **Power analysis reveals that the study is severely underpowered** (power = 0.34–0.68), meaning the experiment has insufficient sample size to reliably detect a true treatment effect if one exists. A non-significant result in an underpowered study is inconclusive, not definitive evidence of no effect.

**Key statistical lessons from this project:**

1. The largest source of variability is between rats (σ²_Rat = 36.06), which dominates Liver-level (σ²_Liver = 14.17) and residual (σ² = 21.17) variation. This confirms that Rat is an essential random effect.
2. Correct nested model specification — explicitly nesting Rat within Treatment — is critical. Both alternative (wrong) models produce biased variance estimates, singular fits, or invalid treatment comparisons.
3. REML should be used for final estimation; ML (REML = FALSE) is required only for likelihood ratio tests involving fixed effects.
4. To increase power and obtain reliable conclusions about treatment effects, future experiments should include more rats per treatment group.

---

## How to Run

1. Clone or download this repository.
2. Ensure `Rats2.xlsx` is in your working directory.
3. Install required R packages (see [Requirements](#requirements)).
4. Open `R12_Nested_Mixed_Model_Analysis_-_Rats_Glycogen_Experiment.ipynb` in **Jupyter Notebook with an R kernel** (e.g., via `IRkernel`) or copy the R code cells into **RStudio**.
5. Run all cells sequentially from top to bottom.

---


