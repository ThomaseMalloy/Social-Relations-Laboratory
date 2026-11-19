# 00_population_values.R
# Population values for simulating multi-informant, multi-round robin (MIRR) family data.
# Sources (Malloy, Pery, & Cook, family views study; positivity and negativity, 1-7 scale):
#   * Variances (actor, partner, relationship, error) and reciprocities: TripleR pooled output
#     across families, FamilyViews_idioUniVcov_29Sep2024 (same values as Table 2 of the paper).
#   * Role means: Table 3 of the paper (each key person's mean actor and partner effects of
#     mother, father, older, younger), reproduced from FamilyViews_wide_idiographic_AP_Effects.
#   * Agreement across informants: Table 4 of the paper (actor and partner effects) and the
#     relationship effects file FAMILY AS UNIT RELATIONSHIP EFFECTS.sav (relationship effects).
#   * Positivity-negativity correlations: actor and partner effects file (-.49 to -.58) and
#     relationship effects file (-.31 to -.42).
# ASSUMED (not available in the files): construct means 5.0 (positivity) and 2.8 (negativity);
#   generalized reciprocity for positivity .90 (TripleR estimate was at the bound, 1.00);
#   both indicators of a construct measure it equally well (error variance from TripleR).
ROLES <- c("M", "F", "O", "Y")                 # mother, father, older child, younger child
ROLE_LABELS <- c(M = "Mother", F = "Father", O = "Older child", Y = "Younger child")
CONSTRUCTS <- c("POS", "NEG")

# ---- variance components (unstandardized, latent; TripleR pooled output) --------------------
pop_var <- data.frame(
  source    = rep(ROLES, each = 2), construct = rep(CONSTRUCTS, 4),
  varA = c(.286423, .270835,  .274565, .185683,  .348410, .337209,  .235311, .266238),
  varP = c(.163149, .099780,  .063968, .045257,  .179150, .140599,  .125617, .103073),
  varR = c(.215526, .267101,  .224228, .329146,  .264017, .346977,  .263796, .270361),
  varE = c(.485524, .364089,  .448454, .325916,  .585744, .434171,  .621973, .465428),
  # generalized (actor-partner) and dyadic reciprocity correlations
  genrec = c(.90, .835170,  .90, .736881,  .90, .697860,  .90, .603691),
  dyrec  = c(.847303, .470231,  .711415, .617145,  .901471, .499374,  .850879, .612300))
pop_var$genrec_TripleR <- c(.998938, .835170, 1, .736881, 1, .697860, 1, .603691)

# ---- role means (Table 3): rows = key person (source); columns = member M, F, O, Y -----------
tab3 <- function(v) matrix(v, 4, 4, byrow = TRUE, dimnames = list(ROLES, ROLES))
role_actor <- list(
  POS = tab3(c(.426, -.015, -.126, -.284,   .419,  .016, -.124, -.311,
               .428, -.041,  .003, -.390,   .362, -.058, -.177, -.127)),
  NEG = tab3(c(-.218, -.114, .051, .281,   -.231, -.081,  .016, .296,
               -.244, -.173, -.022, .439,  -.256, -.134,  .223, .167)))
role_partner <- list(
  POS = tab3(c(.318, -.195, -.019, -.104,   .223, -.090, -.002, -.131,
               .265, -.277,  .164, -.153,   .307, -.243, -.061, -.003)),
  NEG = tab3(c(-.076,  .075, -.137, .139,   -.015, -.033, -.115, .163,
                .009,  .028, -.191, .154,   -.092,  .023, -.017, .086)))
# centre each row exactly (Table 3 entries are rounded)
for (c in CONSTRUCTS) {
  role_actor[[c]]   <- sweep(role_actor[[c]],   1, rowMeans(role_actor[[c]]))
  role_partner[[c]] <- sweep(role_partner[[c]], 1, rowMeans(role_partner[[c]]))
}

construct_mean <- c(POS = 5.0, NEG = 2.8)       # ASSUMED

# ---- positivity-negativity correlations of the random effects ------------------------------
r_PN <- list(actor = -.55, partner = -.55,       # observed -.49 to -.58 (estimates, attenuated)
             rel_same = -.45, rel_cross = -.20)  # observed rel. effects -.31 to -.42

# ---- agreement targets (estimated effects, averaged over roles and constructs) --------------
agree_target <- c(other_other = .48,   # Table 4, other-other, mean of actor and partner
                  self_other  = .36,   # Table 4, self-other
                  relationship = .27)  # .sav file: POS .30, NEG .24

# ---- calibrated shared-view weights (see 03_check_recovery.R) ---------------------------------
# Chosen so that estimated agreement in large samples matches the targets above:
# other-other .49, self-other .38, relationship .28 (POS .24, NEG .31) in 1,500 families.
W_OTHER <- .86   # key person's view of another member's actor/partner effect
W_SELF  <- .02   # key person's view of their own actor/partner effect
W_REL   <- .62   # relationship effects
