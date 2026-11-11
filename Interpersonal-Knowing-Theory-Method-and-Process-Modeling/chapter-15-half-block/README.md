# Chapter 15: Half-block SRM analysis

R code for the half-block design in Chapter 15. In a half-block, actors (rows) rate partners (columns) and partners do not rate actors. Agents rating other agents' introductions, perceivers rating faces, and humans rating digital twins all have this structure. The code is meant to be easier to use than Blocko or general asymmetric block design (ABD) software.

## What it does

- Reads long-format half-block data: one row per rating, with columns `block`, `actor`, `partner`, and the indicators (for example `warm`, `kind`, `capable`, `effective`).
- Estimates actor, partner, and relationship variance in each block with the random-effects ANOVA estimators in Chapter 6 (Equations 6.1 to 6.3), then pools the block estimates (the Blocko approach).
- With two indicators per construct, separates relationship from error. The construct variance of each component is the covariance of that component across the two indicators.
- Estimates actor-actor, partner-partner, and relationship-relationship covariances and correlations between constructs.
- Tests variances (one-tailed) and covariances (two-tailed) with one-sample t tests across blocks.
- Writes actor and partner effects for second-stage modeling.
- Fits the same model in lme4. This handles missing ratings and blocks of different sizes. With two indicators stacked, the dyad term separates relationship from error.
- Analyzes the own-twin design (Table 15.1), in which each human rates their own digital twin and other twins: own-twin vs. other-twin ratings and relationship effects side by side, the own-twin effect and its test, and rater-level values for second-stage modeling.

## Files

| File | What it is |
| --- | --- |
| `half_block_functions.R` | All functions (base R; lme4 only for `hb_lme4()`) |
| `0_make_example_data.R` | Simulates 100 blocks of 4 actors x 4 partners with known values; writes `data/example_half_block.csv` and `data/example_own_twin.csv` |
| `1_half_block_srm.R` | The analysis: pooled estimates, tests, standardized values, correlations, effects, and lme4 |
| `2_own_twin_design.R` | Own-twin effect (Table 15.1 design) |
| `3_verify_estimators.R` | Monte Carlo check of bias and Type I error |
| `output/` | Results from the example data and the verification |

## Steps

Open `chapter-15-half-block.Rproj` in RStudio (or `setwd()` to this folder), then:

```r
source("0_make_example_data.R")   # optional: example data
source("1_half_block_srm.R")      # edit the file name and constructs for your data
source("2_own_twin_design.R")     # own-twin design only
```

To use your own data, change two lines in `1_half_block_srm.R`:

```r
d <- read.csv("your_data.csv")    # columns: block, actor, partner, indicators
constructs <- list(warmth = c("warm", "kind"), competence = c("capable", "effective"))
```

A construct can have one indicator (relationship and error are then combined) or two.

## Notes

- The ANOVA estimators need complete blocks. Use `hb_lme4()` when ratings are missing.
- Blocks may differ in size. Block estimates are averaged without weights, as in Blocko.
- Negative variance estimates are possible when a true variance is near zero. Correlations are not computed when a variance estimate is zero or negative.
- In the own-twin design, the own-twin effect is the relationship effect for the human's own twin minus the mean relationship effect for the other twins. Humans in the same block are not independent, so block means are tested across blocks. In lme4, the own twin is a fixed effect, which keeps it out of the relationship variance.

## Verification

`3_verify_estimators.R` simulated 500 replications of 50 blocks for 4 x 4, 6 x 6, and 3 x 8 blocks. Every variance and covariance estimate was within Monte Carlo error of its true value. With true values of zero, Type I error rates were 4.0% (actor variance), 4.5% (partner covariance), and 5.8% (own-twin effect). See `output/verification.txt`.

## Software

R 4.3.3; lme4 1.1 (optional).
