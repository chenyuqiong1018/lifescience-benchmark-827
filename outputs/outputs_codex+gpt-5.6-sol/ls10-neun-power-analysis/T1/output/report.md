# NeuN power-analysis audit

## Frozen analysis choice

The prompt specifies a two-sided independent t-test, equal group sizes, alpha = 0.05, and power = 0.80. The repeated Sample identifiers were therefore not used to silently substitute a paired design. All NeuN cells are observed; no missing value was imputed.

## Descriptive statistics and effect size

- KD: n = 8, mean = 214.500, sample SD = 10.941402
- CTRL: n = 8, mean = 210.625, sample SD = 22.853180

The equal-variance pooled SD is **17.916224**. In frozen label order (KD minus CTRL), Cohen's d is **0.216284**. The sign records direction; the requested two-sided power calculation uses |d|.

For transparency, pilot-data diagnostics were Shapiro-Wilk p = 0.774894 (KD) and 0.723212 (CTRL), with median-centered Levene p = 0.107338. These small-sample diagnostics do not redefine the test model fixed by the prompt.

## Power result

Direct noncentral-t enumeration finds **337 observations per group**. At 336, power is 0.799373; at 337, it is 0.800542. The statsmodels continuous solution is 336.535394, whose upward rounding independently confirms 337.

## Skill use

`statistical-analysis` supplied the independent-group Cohen's d convention, sample-SD requirement, assumption transparency, and noncentral-t planning model. `code_execution_analysis` was asked for a fixed-data arithmetic cross-check, but its endpoint only returned the code without executing it; the rerunnable local script is therefore the executed source of record. `biomarker_discovery` was opened as required but its external biomarker databases were not queried because they are irrelevant to the supplied measurements and external data are forbidden.
