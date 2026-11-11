# Chapter 14: Triadic idiographic key person design (triadic IKPD)

These scripts do every analysis proposed in Chapter 14. They read a key person's (KP's) triadic assessments of nominees' relationships (Table 14.2), run the social relations model (SRM) in each idiographic round robin with TripleR (or the same equations in base R), and write the grand means, effect estimates, and variance components as data sets for second-stage modeling. Step 2 then tests the hypotheses in Tables 14.6, 14.7, and 14.9.

## Data you need

One round robin per KP per social group (for example friends and coworkers). Each round robin has the KP and at least three nominees. Each cell is the KP's view of how the row person (perceiver) assesses the column person's (target's) relationship with that target's partner. The KP's own row is the KP's own view. The diagonal is not collected.

Either format works:

| Format | Layout |
| --- | --- |
| Long | One row per KP × group × perceiver × target: `kp.id, group, perceiver, target, <construct 1>, <construct 2>, ...` |
| Wide | One row per KP × group, with one column per cell named `<construct>_<perceiver>_<target>`, for example `satisfaction_KP_N1` (the KP's view of N1's satisfaction) and `satisfaction_N1_KP` (the KP's view of how N1 sees the KP's satisfaction). This is how survey software usually exports it. |

The key person must be labelled `KP`. Nominees can have any labels (`N1`, `N2`, `N3` in the example).

## Files

| File | What it does |
| --- | --- |
| `srm_triadic_functions.R` | Base R functions used by all steps: reading and checking the data, the round robin equations (Chapter 4; Warner, Kenny, & Stoto, 1979) for effects, variances, and covariances, and the bivariate covariances. |
| `0_make_example_data.R` | Optional. Simulates 150 KPs, each with 3 friends and 3 coworkers, on listening (X) and relationship satisfaction (Y), from known values. Writes the data in long and wide form and the true values to `output/true_values.txt`. |
| `1_srm_effects.R` | Step 1. Reads your data (set the file and column names at the top), runs the SRM in every round robin, and writes the second-stage data sets below. Uses TripleR when installed and checks it against the base R equations. |
| `2_hypothesis_tests.R` | Step 2. Tests the within-group hypotheses (Table 14.6), the between-group hypotheses (Table 14.7), and the bivariate hypotheses (Table 14.9). Writes `output/triadic_results.txt` and `.csv`. |
| `3_verify_estimators.R` | Optional. Monte Carlo check that the estimators are unbiased for n = 4 and n = 6, and that the tests reject about 5% of the time when the null hypotheses are true. Writes `output/verification.txt`. |
| `4_second_stage_example.R` | Optional. Two second-stage examples with the exported data sets: a model of KP versus nominee target effects, and a lavaan path model of KP egocentric bias across groups. |

## Second-stage data sets (written by Step 1 to `data/`)

| File | One row per | Columns |
| --- | --- | --- |
| `triadic_person_effects.csv` | KP × group × focal person | `kp.id, group, focal, role` (KP or nominee), and for each construct `<c>.a` (triadic perceiver effect α) and `<c>.b` (triadic target effect β) |
| `triadic_relationship_effects.csv` | KP × group × perceiver × target | `kp.id, group, perceiver, target`, and `<c>.g` (relationship effect γ, Table 14.4) |
| `triadic_kp_level.csv` | KP × group | `n`; for each construct: `.mean` (grand mean µ), `.a.kp` and `.b.kp` (KP's perceiver and target effects), `.a.nom` and `.b.nom` (mean nominee effects), `.var.a`, `.var.b`, `.var.g` (perceiver, target, relationship variances), `.var.total` (their sum: total variance, the KP's perceived heterogeneity of the group), `.cov.ab`, `.cov.gg` (perceiver–target and dyadic covariances), `.var.mean` (sampling variance of the grand mean); for each pair of constructs: `.cov.aa`, `.cov.bb`, `.cov.gg.intra`, `.cov.gg.inter` |
| `triadic_kp_wide.csv` | KP | The KP-level file with one column per group, for example `satisfaction.b.kp_friends` and `satisfaction.b.kp_coworkers` |

## Steps

Open `chapter-14-triadic.Rproj` in RStudio, or set the working directory to this folder. Then run:

```r
source("0_make_example_data.R")      # optional: example data with known answers
source("1_srm_effects.R")            # your data: set the file and column names at the top first
source("2_hypothesis_tests.R")       # Tables 14.6, 14.7, 14.9
source("3_verify_estimators.R")      # optional: Monte Carlo check of the estimators and tests
source("4_second_stage_example.R")   # optional: second-stage examples (lme4, lavaan)
```

Packages: TripleR (optional, Step 1), lme4 (Step 2 grand mean test and Step 4), lavaan (Step 4). Everything else is base R.

## How the tests work

- **Unit of analysis.** Each KP gives one estimate of each parameter per group. Tests are t tests across KPs, as for the dyadic IKPD.
- **Perceiver, target, and relationship variances (= 0).** One-sample t tests of the KP variance estimates, one-sided.
- **KP versus average nominee (α_KP = mean α_N, β_KP = mean β_N).** Paired t tests. The egocentric bias hypothesis β_KP > mean β_N is the one-sided version.
- **Between groups.** Paired t tests of grand means, variances, and KP effects across groups, and correlations of the KP effects across groups (consistency).
- **Total variance (perceived heterogeneity).** Each KP's total variance (perceiver + target + relationship, error included with one indicator) measures how much the KP differentiates the relationships in a group. A paired t test (df = K − 1) compares groups. Because each KP is compared with themselves, individual differences in scale use cancel. It is not tested against zero, which it exceeds whenever the ratings are not all identical. The component tests show which source drives a difference.
- **Variance of the grand mean (σ²µ = 0).** A KP has only one grand mean per group, so there is no per-KP variance to put into a t test. Step 2 tests whether KP grand means vary more than sampling error allows: (a) a multilevel model of the raw scores with a KP random intercept and perceivers and targets nested in KPs, likelihood ratio test with the p value halved (the variance is tested at its boundary of zero); (b) a moment estimate, the variance of the KP grand means minus the average sampling variance implied by each KP's SRM components; and (c) a one-sample t test (df = K − 1) of d_k = K/(K − 1)(M_k − mean M)² minus KP k's sampling variance, whose mean is the estimate in (b). A one-sample t test of the grand means themselves tests their mean (for example against the scale midpoint), not their variance.
- **Bivariate (Table 14.9).** The perceiver–perceiver, target–target, and intrapersonal relationship covariances are estimated in each round robin and tested with one-sample t tests across KPs. They are computed with the univariate round robin estimator applied to X + Y and X − Y: Cov = [V(X + Y) − V(X − Y)] / 4, which is exactly the bivariate round robin estimator.

## Testing

- Tested with R 4.3.3 on the example data. `3_verify_estimators.R` shows that every variance and covariance estimator is unbiased for n = 4 and n = 6 and that the tests reject at the nominal 5% rate under the null hypotheses (`output/verification.txt`). `output/triadic_results.txt` recovers the values in `output/true_values.txt`.
- Step 1 checks that the effects sum to zero in every round robin and that grand mean + perceiver + target + relationship effects reproduce every score.
- **TripleR could not be installed in the test environment (CRAN is blocked there), so the TripleR branch of Step 1 has not been run.** With TripleR installed, Step 1 compares TripleR's perceiver and target effects with the base R equations and stops if they differ by more than 0.000001, so any mismatch will show at once.
- With a single indicator, relationship effects include measurement error; with two or more indicators, relationship and error can be separated.
