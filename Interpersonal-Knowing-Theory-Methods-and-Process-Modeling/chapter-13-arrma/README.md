# Chapter 13: Estimating ARRMA from TripleR output

These scripts read actor, partner, and relationship effects from a TripleR analysis of a round robin. They then estimate the minimal individual and dyadic ARRMA models of Chapter 13. Each model is estimated two ways: with the standard equations and with the matrix equations.

| ARRMA parameter | Individual (Eq. 13.5) | Dyadic (Eq. 13.12) |
| --- | --- | --- |
| A, reciprocity | cor(α_r, β_r) | cor(γ_ij, γ_ji) |
| B, unique assumed reciprocity b* | Eq. 13.10 | Eq. 13.17 |
| C, unique meta-accuracy c* | Eq. 13.11 | Eq. 13.18 |
| Matrix form | θ̂ = R_zz⁻¹ r_zy, Z = [α_r β_r], y = α_mp | θ̂ = R_zz⁻¹ r_zy, Z = [γ_ij γ_ji], y = γ_imp,ji |

## Data you need

A long round robin file with one row per directed pair. It needs a round robin id, an actor (perceiver) id, a partner (target) id, the actor's **judgment** of the partner, and the actor's **metaperception**: the actor's prediction of how the partner judged the actor.

## Files

| File | What it does |
| --- | --- |
| `0_make_example_data.R` | Optional. Simulates 40 round robins of 6 people from known ARRMA values. It also computes the SRM effects in base R with the Chapter 4 equations (4.1 to 4.3), so Step 2 can run without TripleR. |
| `1_TripleR_export_effects.R` | Step 1. Runs TripleR on the judgment and the metaperception. It writes the person file and the relationship file used by Step 2. Needs TripleR. |
| `2_arrma_estimate.R` | Step 2. Reads the two effects files and estimates individual and dyadic ARRMA with the standard equations and the matrix equations. It checks that both methods agree with each other and with `lm()`, then reports tests and bootstrap confidence intervals. Base R only. |
| `data/arrma_person_effects.csv` | `group.id, id, judgment.a, judgment.p, meta.a, meta.p` |
| `data/arrma_relationship_effects.csv` | `group.id, actor.id, partner.id, judgment.r, meta.r`. For actor *i* and partner *j*, `meta.r` is γ_imp,ji: *i*'s prediction of *j*'s unique response to *i*. |
| `data/example_round_robins.csv` | Simulated example: 40 round robins of 6 people, 1,200 directed ratings. Columns `group.id, actor.id, partner.id, judgment, meta`. True values: individual b* = .50, c* = .25, reciprocity = .30; dyadic b* = .45, c* = .20, reciprocity = .40. |
| `output/arrma_results.txt` | Results for the example data, for comparison |
| `output/arrma_apa_table.rtf` | APA-style table of the ARRMA estimates with 95% CIs and *p* values. It opens in Word. A plain text version (`.txt`) and a `.csv` are written too. |

## Steps

Open `chapter-13-arrma.Rproj` in RStudio, or set the working directory to this folder, for example with `setwd("C:/.../chapter-13-arrma")`. Then run:

```r
source("0_make_example_data.R")       # optional: example data with known answers
source("1_TripleR_export_effects.R")  # your data: set the column names at the top first
source("2_arrma_estimate.R")          # ARRMA, standard and matrix equations
```

Already have TripleR output? Skip Step 1. Save the person effects (`$effects`) and the relationship effects (`$effectsRel`) for the judgment and the metaperception as CSV files. Then type your file and column names into the settings block at the top of `2_arrma_estimate.R`.

## Notes

- The standard equations and the matrix solution θ̂ = R_zz⁻¹ r_zy give identical b* and c*. The script stops with an error if they differ.
- **Individual level:** t tests use df = N − G − 2. SRM effects sum to zero within each round robin, which costs one degree of freedom per round robin, and the model has 2 predictors.
- **APA table:** CIs are percentile bootstrap intervals from resampling round robins. Individual-level *p* values are t tests (df = N − G − 1 for correlations, N − G − 2 for b* and c*). Dyadic *p* values are two-sided bootstrap *p* values.
- **Dyadic level:** each dyad enters twice, as (*i*, *j*) and as (*j*, *i*), and people appear in several dyads. Use the bootstrap confidence intervals, which resample whole round robins, rather than `lm()` standard errors.
- **Single indicators:** relationship effects from one indicator contain error, which pulls dyadic correlations toward zero. With two or more indicators, relationship and error can be separated. In the example data, the dyadic estimates are a little below the true values for this reason.
- **Testing:** the scripts were tested with R 4.3.3 on the example data. TripleR could not be installed in the test environment, so Step 1 was not run there. It uses the same `RR()` call and `$effects` layout as the Chapter 10 scripts.
