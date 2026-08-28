# Treatment-response logistic model

## Coding and analysis set

Treatment response was the modeled event: `Efficacy = PR` was coded 1, and `SD` or `PD` was coded 0. Female was the Gender reference level; `Gender_Male` equals 1 for Male. Complete cases were defined only over Efficacy, BMI, Age, and Gender. All 80 candidate rows were complete, leaving n = 80.

The prespecified multiple logistic regression was `logit(P(PR)) = intercept + BMI + Age + Gender_Male`. No distractor covariate from the workbook entered the model. Maximum likelihood converged. The overall likelihood-ratio test versus the intercept-only model was chi-square(3) = 13.696324, p = 0.00334905.

## Age coefficient

Holding BMI and gender constant, the **Age log-odds coefficient was -0.079508** (SE = 0.026298, Wald z = -3.023397, **two-sided p = 0.00249954**). The corresponding odds ratio was 0.923570 per age unit (Wald 95% CI 0.877173 to 0.972422).

`statistical-analysis` guided selection of multiple logistic regression for a binary outcome, explicit outcome/reference coding, complete-case disclosure, coefficient/SE/z/p/OR reporting, and convergence/model-fit checks.
