# Treatment-response logistic regression

## Model definition

The binary modeled event is treatment response: **Efficacy = PR was coded 1**, while **SD and PD were coded 0**. Gender used **Female as the reference level** and `Gender_Male = 1` for Male. The logit model included only BMI, Age, and Gender_Male, plus an intercept; distractor workbook covariates were excluded.

Complete-case filtering was restricted to Efficacy, BMI, Age, and Gender. All 80 data rows were complete, so the fitted sample was n = 80 with 41 modeled responses and 39 non-responses.

## Requested age result

Holding BMI and gender fixed, the **Age log-odds coefficient is -0.079508** (SE = 0.026298, z = -3.023397, **two-sided p = 0.00249954**). Its odds ratio is 0.923570 per one-unit increase in age.

The model converged using maximum-likelihood binary logistic regression with a logit link. Coefficients, standard errors, Wald z statistics, two-sided p-values, and exponentiated coefficients are provided in `model_coefficients.csv`.
