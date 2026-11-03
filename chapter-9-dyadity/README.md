# Dyadity analysis in R (Chapter 9)

`dyadity_analysis.R` runs every analysis in Chapter 9 in base R; no packages are needed.

1. **SRM effects** for each round robin (Eqs. 9.3–9.5).
2. **Expressed and elicited dyadity** for every member (Eqs. 9.6–9.9). With two or more indicators of the same construct, these use cross-products averaged over pairs of indicators (Eqs. 9.6a and 9.8a).
3. **Expressed minus elicited dyadity** for each member, ∆D (Eq. 9.10).
4. **Monte Carlo tests** of individual differences in expressed and elicited dyadity.
5. **Role-specific tests** of ∆D: a one-sample t test for one role, or a paired t test for two roles.
6. **The key-person contrast** for egocentric bias in the IKPD (Eqs. 9.11 and 9.12).
7. **Matrix-algebra verification.** `verify_matrix()` recomputes Eqs. 9.3–9.10 with the matrix formulas and checks the identities for each round robin: zero row and column sums of Γ̂, D̄_Act = D̄_Prt, 1′∆d = 0, and C·C⁻¹ = I.

## Steps

```r
source("dyadity_analysis.R")
run_example()      # reproduces the Chapter 8 round robins; takes a few seconds

res <- dyadity_analysis(mydata, group = "family", actor = "rater",
                        partner = "target",
                        y = "listen1",                  # one indicator (univariate)
                        # y = c("listen1", "listen2"),  # two or more indicators (multivariate)
                        B = 9999, seed = 1)
print(res)                                   # person table, group table, Monte Carlo tests
verify_matrix(res)                           # matrix-algebra checks
write.csv(res$persons, "dyadity_by_person.csv", row.names = FALSE)

roles <- data.frame(family = ..., person = ..., role = ...)   # names: group, person, role
role_test(res, roles, role = "mother")                       # Delta D for mothers vs 0
role_test(res, roles, role = "mother", role2 = "father")     # mothers vs fathers (paired)
kp_contrast_test(res, roles, kp_role = "KP")                 # IKPD egocentric bias
```

## Data format

The data are in long format, with one row per directed dyad. Each row holds:

- a group identifier
- the actor (the person who rates or behaves)
- the partner
- the score (with two or more indicators of the same construct, one column per indicator)

Each group must be a complete round robin with at least four members. Groups can differ in size.

The `roles` data frame needs three columns named `group`, `person`, and `role`. Its `group` and `person` values must match those in the data.

## Reading the output

- **D_exp and D_elic** are each member's expressed and elicited dyadity. `ratio_exp` and `ratio_elic` are the same indices divided by the round-robin mean, so 1.00 means average.
- **Delta_D** is expressed minus elicited dyadity. It always sums to 0 within a round robin, so only its variance, or its value for a specific role, can be tested.
- **The Monte Carlo test** compares the mean within-group variance of the indices (`observed_T`) with its expected value when nobody differs (`null_mean`). It does not compare with 0. The p value is (1 + number of simulated values at least as large as the observed value) / (B + 1).

## Practical notes

- **The Monte Carlo test is calibrated.** With simulated data in which nobody truly differs, it rejected at α = .05 in 3–5% of samples with one indicator and 2–8% with two or three indicators.

- **Power depends on round-robin size.** In one simulation, a member had four times the relationship variance of the others. The test detected it 82% of the time with 20 eight-person round robins, but only 13% of the time with 50 four-person round robins.
- **Two or more indicators remove measurement error** from the dyadity indices. A member's index can then be negative in a sample, which is read as no dyadity beyond error.
- **With one indicator and n = 4,** the expressed and elicited variance tests give identical results.
