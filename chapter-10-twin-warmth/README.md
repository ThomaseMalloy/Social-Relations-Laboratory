# Chapter 10: Human agreeableness and digital twin warmth

This folder reproduces the digital twin example in Chapter 10: the first-stage SRM analysis (Table 10.1) and the multiple-group SEM (Table 10.2, Table 10.3, Figure 10.6). The data are from Malloy and Dumais (2026), distinct-persona condition. Each twin was instantiated from a different real person in the Twin-2K-500 data set (Toubia et al., 2025).

## Design

| | Kimi K2.6 | DeepSeek V4 Pro | Nemotron 3 Ultra |
| --- | --- | --- | --- |
| Round robins (n = 5) | 50 | 50 | 50 |
| Twins | 250 | 250 | 250 |
| Directed ratings | 1,000 | 1,000 | 1,000 |

After each interaction a twin rated its partner on **warm** and **kind** (1 = low, 5 = high).

## Files

| File | What it is |
| --- | --- |
| `data/kimi_distinct_round_robins.csv` | Raw round robin ratings, Kimi. One row per directed pair. Key columns: `roundrobin`, `actor.id`, `partner.id`, `warm`, `kind`. The other columns record the model, run, and seed. |
| `data/deepseek_distinct_round_robins.csv` | Same for DeepSeek |
| `data/nemotron_distinct_round_robins.csv` | Same for Nemotron |
| `data/human_agreeableness.csv` | Agreeableness (Big Five, 1 to 5) of the real person behind each twin: `platform`, `group.id`, `id`, `agreeableness` |
| `data/ch10_twin_warmth_srm_effects.csv` | First-stage actor and partner effects for warm and kind, one row per twin, with agreeableness. This is the lavaan input. |
| `1_TripleR_srm_effects.R` | Step 1. TripleR SRM analysis of warm and kind on each platform; writes the effects file |
| `2_lavaan_sem.R` | Step 2. Multiple-group SEM (Models 1 to 5), difference tests, and final estimates |
| `3_variance_components_table_10_1.R` | Table 10.1 variance partitioning in base R, plus a check that the effects file matches the raw ratings |
| `output/` | Output of scripts 2 and 3 for comparison |

## Steps

Open `chapter-10-twin-warmth.Rproj` in RStudio (this sets the working directory), or set it yourself with `setwd("path/to/chapter-10-twin-warmth")`. Then run:

```r
source("1_TripleR_srm_effects.R")              # needs TripleR; rewrites data/ch10_twin_warmth_srm_effects.csv
source("2_lavaan_sem.R")                       # needs lavaan; Tables 10.2 and 10.3, Figure 10.6
source("3_variance_components_table_10_1.R")   # base R; Table 10.1
```

Step 1 is optional because its output, the effects file, is included. Script 3 checks that file against the raw ratings without TripleR.

## Model

```
WA =~ warm.actor + kind.actor            # latent twin warmth actor effect
WP =~ warm.partner + kind.partner        # latent twin warmth partner effect
warm.actor ~~ warm.partner               # same-item residual covariance
WA ~ agreeableness
WP ~ WA + agreeableness                  # direct path = accuracy
```

Platform is the group. Standard errors are cluster-robust by round robin, and nested models are compared with the scaled chi-square difference test (Satorra & Bentler, 2001). In the final model (Model 5), the accuracy path and the residual covariance are equal across platforms: χ²(10) = 4.55, p = .92, RMSEA = .00, CFI = 1.00, accuracy b = .27 (SE = .03).

lavaan may warn that the covariance matrix of the parameter estimates is not positive definite. This follows from the effects summing to zero within each round robin, which the cluster-robust standard errors take into account. The estimates are not affected.

## Software

R 4.3.3, lavaan 0.6-17, TripleR (Schönbrodt, Back, & Schmukle, 2012).

## References

Malloy, T. E., & Dumais, T. (2026). *Digital twins' dyadic interactions: Convergence of inter-persona assessments and human origins* [Unpublished manuscript]. Rhode Island College.

Rosseel, Y. (2012). lavaan: An R package for structural equation modeling. *Journal of Statistical Software, 48*(2), 1–36. https://doi.org/10.18637/jss.v048.i02

Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test statistic for moment structure analysis. *Psychometrika, 66*(4), 507–514. https://doi.org/10.1007/BF02296192

Schönbrodt, F. D., Back, M. D., & Schmukle, S. C. (2012). TripleR: An R package for social relations analyses based on round-robin designs. *Behavior Research Methods, 44*(2), 455–470. https://doi.org/10.3758/s13428-011-0150-4

Toubia, O., Gui, G. Z., Peng, T., Merlau, D. J., Li, A., & Chen, H. (2025). Database report: Twin-2K-500: A data set for building digital twins of over 2,000 people based on their answers to over 500 questions. *Marketing Science, 44*(6), 1446–1455. https://doi.org/10.1287/mksc.2025.0262

Warner, R. M., Kenny, D. A., & Stoto, M. (1979). A new round robin analysis of variance for social interaction data. *Journal of Personality and Social Psychology, 37*(10), 1742–1757. https://doi.org/10.1037/0022-3514.37.10.1742
