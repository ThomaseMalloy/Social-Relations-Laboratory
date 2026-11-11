# Chapter 8: Multi-informant, multi-round robin (MIRR) families — simulation and TripleR

R code to simulate MIRR data for four-person families (mother, father, older child, younger
child), analyze every key person's idiographic round robin with TripleR, and save all
first-stage SRM parameters in long and wide formats for second-stage analyses in R, SPSS,
or other software.

**The data are SIMULATED.** Population values come from the family views study
(Malloy, Pery, & Cook): positivity and negativity, two indicators each, 1–7 scale.

## Design

Each family member is a key person (source) who rates all 12 directed relationships among the
four members: their own behavior toward others (self-to-other), others' behavior toward them
(other-to-self), and others' behavior toward one another (triadic). Each family therefore
yields four idiographic 4 × 4 round robins. With 50 families there are 50 round robins per
role and 200 in all.

## Files

| File | Purpose |
|---|---|
| `00_population_values.R` | Population values and their sources (see below). |
| `sim_functions.R` | Generating model. |
| `srm_functions.R` | Built-in round-robin SRM estimator following TripleR's formulas (used when TripleR is not installed, and to check TripleR). |
| `param_functions.R` | Variable naming, long and wide files, CSV/SPSS/R output, agreement. |
| `01_simulate.R` | Simulates 50 families → `data/mirr_raw_long` and `data/mirr_raw_wide`. |
| `02_srm_analysis.R` | TripleR for every family × key person (latent constructs, plus positivity × negativity), prints the output, and writes the parameter files. Also runs TripleR's pooled multi-group analysis for each key-person role. |
| `03_second_stage.R` | Second-stage examples: Table 2 (variance components, t tests across families), Table 3 (role means of actor and partner effects), Table 4 (self–other and other–other agreement). |
| `04_check_recovery.R` | Simulates 2,000 families and checks that the population values are recovered. |
| `05_compare_engines.R` | With TripleR installed: compares every TripleR estimate with the built-in estimator. |

Run in order: `Rscript 01_simulate.R`, `Rscript 02_srm_analysis.R`, `Rscript 03_second_stage.R`.
Packages: TripleR (`install.packages("TripleR")`), haven (for SPSS files), Matrix, MASS.

## Output files

| File | Contents |
|---|---|
| `data/mirr_raw_long.*` | One row per family × source × actor × partner: POS1, POS2, NEG1, NEG2. This is the layout TripleR reads. |
| `data/mirr_raw_wide.*` | One row per family; 192 ratings named `SOURCE_ITEM_ACTORPARTNER`. |
| `output/tripler_output_by_family.txt` | Printed TripleR output for all 200 round robins. |
| `output/tripler_output_pooled.txt` | Printed TripleR multi-group output, one per key-person role and construct. |
| `output/srm_params_long.*` | One row per family × source × parameter (× target): family, source, variable, target, estimate, label, name. |
| `output/srm_params_wide.*` | One row per family: all 312 parameters, mother's block, then father's, then older child's, then younger child's. The order is identical for every family. |
| `output/variable_dictionary.csv` | Every wide variable name with its meaning (also the SPSS variable labels). |
| `output/table2_*.csv`, `table3_*.csv`, `table4_*.csv` | Second-stage results. |

Each file is written as `.csv` (any software), `.sav` (SPSS, with variable labels), and
`.rds` (R).

## Variable names: SOURCE_VARIABLE_TARGET

SOURCE is the key person who made the ratings: M, F, O, Y. VARIABLE is construct + parameter.
TARGET is the member (M, F, O, Y) or directed dyad (FO = father toward older child). Names
for parameters of the whole round robin have no target.

| Example | Meaning |
|---|---|
| `M_POS1_FO` | Mother's rating, positivity item 1, father toward older child (raw data) |
| `M_POSmean` | Mother's within-person grand mean, positivity |
| `M_POSact_F` / `M_POSprt_F` | Mother's view of father's actor / partner effect |
| `M_POSrel_FO` | Mother's view of the father-to-older-child relationship effect |
| `M_POSvarA`, `varP`, `varR`, `varE` | Actor, partner, relationship, error variance |
| `M_POSstdA` … `stdE` | Variance as a proportion of total |
| `M_POScovAP`, `M_POSgenrec` | Actor–partner covariance and generalized reciprocity r |
| `M_POScovR`, `M_POSdyrec` | Relationship covariance and dyadic reciprocity r |
| `M_PNcovAA`, `M_PNrAA` … | Positivity × negativity covariances and correlations (actor–actor, partner–partner, actor–partner, partner–actor, relationship same and reversed direction) |

Within each source block the order is: POS, then NEG (mean, variances, covariances,
standardized, reciprocities, actor effects M F O Y, partner effects M F O Y, relationship
effects MF MO MY FM FO FY OM OF OY YM YF YO), then the positivity × negativity block.

## Population values

| Quantity | Source |
|---|---|
| Actor, partner, relationship, error variances | TripleR pooled output (FamilyViews_idioUniVcov, = Table 2) |
| Generalized and dyadic reciprocity | Same file. Generalized reciprocity for positivity was at the bound (1.00) for every role; **set to .90** |
| Role means (each key person's view of each member's actor and partner effect) | Table 3 (reproduced from the actor–partner effects file) |
| Agreement between informants | Table 4 (actor, partner) and the relationship effects file (.30 positivity, .24 negativity) |
| Positivity–negativity correlations | Actor and partner effects file (−.49 to −.58); relationship effects file (−.31 to −.42) |
| Construct means | **Assumed**: positivity 5.0, negativity 2.8 (no raw means were available) |
| Indicators | Both indicators measure the construct equally well (**assumed**); each has its own error |

## Generating model

Rating by key person k of member i's behavior toward member j:

X = mean + role actor mean_k(i) + role partner mean_k(j) + a_k(i) + b_k(j) + g_k(i→j) + error

Agreement between informants comes from a family-level shared effect that every key person
sees, mixed with each key person's own view: effect_k = √w · shared + √(1 − w) · own_k.
The weights (w = .86 for views of other members, .02 for views of oneself, .62 for
relationship effects) were calibrated so that the estimated agreement in large samples matches
Table 4 and the relationship-effects file. The mix leaves each key person's variances,
reciprocities, and positivity–negativity correlations at their population values. Ratings are
rounded to integers and kept within 1–7; the error variance is reduced by 1/12 before rounding
(Sheppard's correction), so the rounded data reproduce the error variance.

## Notes for users

* Negative variances and correlations outside ±1 can occur and are kept
  as estimated (correlations are left empty when a variance is not positive). Inferencesnces are
  at the second stage, across families (`03_second_stage.R`).
* The TripleR route in
  `02_srm_analysis.R` is written against TripleR's documented output (`varComp`, `effects`,
  `effectsRel`, `bivariate`) but has not been run. The files supplied were produced with the
  built-in estimator. Run `02_srm_analysis.R` and then `05_compare_engines.R` with TripleR
  installed to confirm the two agree.
* `04_check_recovery.R` results (2,000 families) are in `output/recovery_check.csv`.
