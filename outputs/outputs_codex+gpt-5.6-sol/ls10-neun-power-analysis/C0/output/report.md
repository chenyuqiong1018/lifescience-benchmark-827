# NeuN effect-size and power analysis

The two frozen groups are **KD** and **CTRL**, with 8 and 8 non-missing observations. Their means are 214.500 and 210.625; their sample standard deviations (`ddof=1`) are 10.941402 and 22.853180.

Using the usual equal-variance pooled standard deviation,

`s_p = sqrt(((n1-1)s1^2 + (n2-1)s2^2) / (n1+n2-2)) = 17.916224`.

With the frozen group order KD minus CTRL, Cohen's d is **0.216284**. The positive sign means the observed KD mean is higher; power uses `abs(d)` because the requested test is two-sided.

For an equal-size, two-sided independent t-test with alpha 0.05 and target power 0.80, the continuous solution is 336.535394 observations per group. Rounding upward gives **337 observations per group**. The resulting power is 0.800542; 336 per group would give 0.799373, so 337 is the smallest integer meeting the target under this model.
