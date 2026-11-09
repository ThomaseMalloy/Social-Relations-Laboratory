# Chapter 15: Asymmetric block design (ABD) SRM analysis

R code for the asymmetric block design in Chapter 15. In an ABD block there are two groups (for example, 3 high-status and 3 low-status agents). Members of each group rate every member of the other group, and no one rates members of their own group. The code reproduces Blocko (Kenny & Xuan) in base R, so readers do not need Fortran.

## What it does

- Reads long-format ABD data: one row per directed rating, with columns `block`, `actor`, `actor_group`, `partner`, and the indicators (for example `warm`, `kind`, `capable`, `effective`).
- For each direction (for example High→Low and Low→High), estimates actor, partner, and relationship variance in each block with the random-effects ANOVA estimators (Chapter 6, Eqs. 6.1 to 6.3), then pools the block estimates across blocks.
- With two indicators per construct, separates stable construct variance (the covariance across indicators) from unstable variance, and reports percent stable variance.
- Estimates reciprocity: the actor–partner correlation for each group (does how a person rates the other group go with how the other group rates them?) and the relationship (dyadic) correlation (does i's unique rating of j go with j's unique rating of i?).
- Tests variances (one-tailed) and covariances (two-tailed) with one-sample t tests across blocks.
- Reports relative variances as Blocko does: the share of the total, with negative estimates set to zero.
- Writes actor and partner effects for second-stage modeling, and compares conditions measured on the same blocks (for example, 0 vs. 2 conversational turns) with paired t tests (`abd_compare()`).

## Files

| File | What it is |
| --- | --- |
| `abd_functions.R` | All functions (base R) |
| `data/claude_status_abd.csv` | The Chapter 15 social-status study: 50 blocks of 3 high- and 3 low-status agents created by Claude, 900 directed ratings at zero acquaintance |
| `1_abd_analysis.R` | The analysis: variances, stable variance, reciprocity, tests, effects |
| `2_compare_with_blocko.R` | Checks the results against Blocko's output for the same data |
| `3_verify_estimators.R` | Monte Carlo check of bias and Type I error |
| `output/` | Results, the Blocko comparison, and the verification |

## Steps

Open `chapter-15-abd.Rproj` in RStudio (or `setwd()` to this folder), then:

```r
source("1_abd_analysis.R")         # the analysis
source("2_compare_with_blocko.R")  # optional: check against Blocko
```

To use your own data, change three lines in `1_abd_analysis.R`:

```r
d <- read.csv("your_data.csv")    # columns: block, actor, actor_group, partner, indicators
constructs <- list(warmth = c("warm", "kind"), competence = c("capable", "effective"))
groups <- c("High", "Low")        # the two values of actor_group
```

Output terms are named `direction.variable.component`. For example, `High_Low.warmth.partner` is the stable partner variance of warmth when high-status agents rate low-status agents, and `recip.warmth.High` is the actor–partner reciprocity of high-status agents.

## Agreement with Blocko

On the Chapter 15 data, all 70 values compared (absolute variances for each indicator and construct, relative stable construct variances, percent stable variance, and all 18 reciprocity correlations) are identical to Blocko's output to three decimals. See `output/comparison_with_blocko.txt`.

## Verification

`3_verify_estimators.R` simulated 300 replications of 50 blocks for blocks of 3 + 3 and 4 + 6. Every variance and reciprocity covariance was within Monte Carlo error of its true value. With true values of zero, Type I error rates were 3.4% (actor variance), 4.1% (actor–partner covariance), and 4.1% (relationship covariance). See `output/verification.txt`.

## Notes

- Blocks must be complete, as in Blocko. Groups may differ in size (for example, 4 + 6). For missing ratings, analyze each direction with `hb_lme4()` in the Chapter 15 half-block folder.
- Each direction of an ABD is a half-block. The variance estimates for one direction are the same as those from the half-block code; the ABD code adds the reciprocity between directions.

## Software

R 4.3.3 (base R only).
