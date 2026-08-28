# NeuN standardized difference and power

## Data handling

The frozen file contains two labels in first-occurrence order: **KD** and **CTRL**. Missing NeuN cells would be omitted rather than converted to zero; none are missing, leaving n = 8 and n = 8. Although sample identifiers repeat across labels, the requested design is explicitly a two-sided **independent** t-test, so no paired analysis was substituted.

## Effect size

KD: M = 214.500, sample SD = 10.941402.  
CTRL: M = 210.625, sample SD = 22.853180.

The equal-variance pooled SD is **17.916224**. Cohen's d, defined in the frozen label order as `(KD mean - CTRL mean) / pooled SD`, is **0.216284**. Its sign therefore records direction; the two-sided power calculation uses its absolute magnitude.

## Required sample size

Using the noncentral-t independent-samples model with equal group sizes, alpha = 0.05, target power = 0.80, and a two-sided alternative, the continuous solution is 336.535394 per group. Rounding upward yields **337 per group**. Power is 0.800542 at 337, while 336 gives 0.799373; therefore 337 is the minimum integer meeting 0.80 under the specified model.

The `statistical-analysis` skill guided the independent-group Cohen's d definition, use of sample SDs, statsmodels power solver, transparent direction convention, and upward/minimality check. No inferential group-comparison p-value was requested or added.
