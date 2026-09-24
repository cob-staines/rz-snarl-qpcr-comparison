# RZ Lab vs. SNARL Bd qPCR comparison README.md

## Objective
Analyze Bd qPCR results from a crossed experimental design to identify quantitative biases between extraction methods and qPCR methods between the two labs (RZ & SNARL), ideally to arrive at a formula to make quantities from the two labs comparable (if they are not already).

## Data

### Experiment data
`qpcr_compare_data_260904.xlsx` (in `$data_dir/bd_qpcr_results/rz_snarl_qpcr_comparison/`, not tracked in git). This file supersedes the earlier RZ plate results file (`RIBBiTR_SierraNevadaComparision_..._rerun2_Diluted_results.xlsx`), which is no longer used.

82 frogs (61 from Pitt/RZ lab, 21 from Sierra nevada/SNARL) were swabbed in duplicate (i.e. 2 swabs, A & B) by the collection lab. Swab A had DNA extraction done at collection lab. The extract was then split in 2, half was run though qPCR at the collection lab. Swab B and the other half of Swab A extract were sent to the other lab, where swab B was extracted and both extracts were run through qPCR. This results in 3 qPCR results for each frog (except in exceptions where a swab or extract was lost/not run).

| Result | Swab | Extraction (lab) | qPCR lab |
|---|---|---|---|
| 1 | A | collection lab | collection lab |
| 2 | A (same extract as 1) | collection lab | other lab |
| 3 | B | other lab | other lab |

- 1 vs 2 isolates the **qPCR** effect (identical extract).
- 2 vs 3 isolates the **extraction** effect, confounded with swab-to-swab variability.

Because of the crossed design, each pairwise comparison draws on a different set of frogs (collection lab = population):

| Panel | Comparison | Frogs | Shared extract? |
|---|---|---|---|
| A | Qiagen-RZ vs Qiagen-SNARL | RZ-collected (Pittsburgh) | yes (qPCR effect) |
| B | PrepMan-RZ vs PrepMan-SNARL | SNARL-collected (Sierra Nevada) | yes (qPCR effect) |
| C | Qiagen-RZ vs PrepMan-RZ | SNARL-collected | no |
| D | Qiagen-SNARL vs PrepMan-SNARL | RZ-collected | no |
| E | Qiagen-RZ vs PrepMan-SNARL | both | no |

Lab, population, and load level are therefore confounded when comparing panels (e.g. C vs D).

Notes:
- Extraction method maps 1-to-1 to lab: **RZ = Qiagen DNeasy**, **SNARL = PrepMan**.
- Both labs quantify against ITS1 plasmid standards, but from different standard sources.
- **`bd_load` (ITS1 copies per swab) is the modeled quantity.** It is corrected for dilution, which differed between the two labs; raw per-reaction quantities (`bd_start_quant`) are not directly comparable.
- No technical (well) replicates are available in the experiment data.
- Many results are zero. A zero means "not detected", which may reflect noise or low detection sensitivity, not necessarily true absence. Non-detects have `bd_start_quant = NA` and `bd_load = 0`.

### Replicate data
Pulled from the RIBBiTR database (`survey_data.bd_qpcr_results`, joined to `sample` and `capture`), predating the experiment, to quantify within-lab variability:
- **qPCR replicates**: the same extract (`sample_name_bd` + `extraction_plate_name`) run more than once within one qPCR lab.
- **Swab replicates**: multiple swabs from the same capture (`capture_id`) through the same extraction and qPCR lab, using one result per swab.

Only ITS1 results from the two labs that pass IPC (or have no IPC) are used. Zeros are explicit in `bd_its1_copies_per_swab`, so results where it is NA are dropped. Experiment swabs are excluded so they are not counted twice.

Replicates may be pooled across species, sites, and time. Pooling across labs is under consideration, pending visual exploration.

Findings so far:
- **Swab replicates exist only for RZ** (Qiagen extraction, RZ qPCR; source projects/populations to be confirmed); there are none for SNARL. Within-lab swab variability cannot be compared between labs, and swab noise for the Sierra population is only indirectly available from the experiment (panel C).
- qPCR replicates are few. Same-plate pairs are likely technical well replicates (see `replicates` in `bd_qpcr_results`) rather than re-runs.
- RZ swab replicates initially looked more variable than the experiment, but their spread among both-positive pairs (SD of log10 difference ≈ 0.82) is similar to experiment panel D (≈ 0.79, Pittsburgh population), and higher than panel C (≈ 0.43, Sierra population, ~10 both-positive pairs). This suggests noise depends on load/population rather than a systematic problem with the swab replicates, which are retained for now.
- Difference vs. mean (Bland–Altman) plots show the noise-vs-load curves of the different-swab comparisons (experiment C–E and RZ swab replicates) line up: spread depends on load, not population. **The RZ swab replicates will inform a shared, load-dependent swab-noise term.**

## Hypothesis
qPCR results from the SNARL lab underestimate low Bd quantities compared to the RZ qPCR protocol, but are otherwise comparable.

## Methods
Bayesian hurdle-lognormal models in brms (cmdstanr backend) on Bd load (ITS1 copies per swab). Each unit (replicate group, or frog) has a latent load expressed as a group-level intercept in `mu`, correlated with a group-level intercept in `hu`; the correlation captures load-dependent detection (zeros treated as non-detections).

Bd load within a population is generally modeled as a hurdle-lognormal quantity. Here, the goal is not to represent a population, but to compare relative values between control groups. The model must therefore handle zeros as non-detections (possibly load-dependent detection), rather than as true absence.

Two-stage approach:
1. **Model 1: Replicates**. Correlated group-level intercepts per replicate group in `mu` and `hu`; qPCR lab effect on `hu`; `sigma` by noise type (swab replicates: total noise; qPCR replicates: noise between qPCR re-runs). Posterior draws and summaries are exported for use as priors in Model 2.
   - Kept simple given limited replicates: a first version with a qPCR run (extract × plate) effect and a `sigma` group-level intercept (load-dependent noise) had 21 divergent transitions; neither term was identified (nearly every qPCR run was a single well; 2–3 positives per group).
   - First-version results: detection increases with load (`cor(mu, hu)` ≈ −0.55, 95% CI −0.91 to −0.04); swab noise ≈ 0.50 log10 per result; qPCR noise and SNARL detection (`hu_qpcr_labsnarl` −0.81 ± 1.84) are poorly informed by replicates and must be learned mainly from the experiment. Replicates cannot answer Key question 1.
   - Group-level detection estimates fall on distinct strands by detected fraction (e.g. 2/2, 1/2, 0/2): with 2–3 results per group, the detected fraction dominates each group's `hu` estimate, and the `mu`–`hu` correlation adds a load-dependent slope within each strand. All-zero groups cluster together, with latent loads informed only by the prior and the correlation. The population-level detection curve implied by the model lies between the strands and is the quantity relevant to Model 2.
2. **Model 2: Experiment**. brms non-linear hurdle-lognormal. For positive results, expected log load `mu = L + d + s · SNARL_qPCR · (L − L_ref)`:
   - `L`: latent swab load = intercept + frog-level intercept (true load) + swab-level intercept (swab + extraction noise, shared by results 1 & 2 which share an extract; SD pooled across extraction methods). Separate SDs per extraction method were strongly negatively correlated in the posterior: each frog has exactly one Qiagen and one PrepMan swab, so only their combined variance is identified. **The design cannot tell whether one extraction method is noisier than the other** (variability part of Key question 3).
   - `d`: constant biases for PrepMan extraction (vs Qiagen) and SNARL qPCR (vs RZ); reference pipeline is Qiagen + RZ qPCR.
   - `s`: SNARL qPCR slope deviation; SNARL − RZ = `d_snarl + s·(L − L_ref)`. `s > 0` = SNARL increasingly underestimates below `L_ref` (the hypothesis); `s = 0` = constant bias. `L_ref` = median positive load in the experiment.
   - `hu` (detection): intercepts by extraction method and qPCR lab, plus a frog-level intercept correlated with the frog-level intercept in `L` (correlation approach, as in Model 1). Labs differ in baseline detection but share the detection-vs-load slope; explicit lab-specific detection curves were not used.
   - `sigma`: residual (qPCR) noise, pooled across labs. A first version with `sigma` by qPCR lab (SNARL ≈ 0.38, RZ ≈ 0.16 log10) was weakly identified (the design mainly identifies the sum of the two labs' qPCR noise; RZ `sigma` ESS 205; 5 divergences), so lab-specific qPCR noise is exploratory only.
   - **No population (`collect_lab`) effects**: the focus is on extraction and qPCR lab differences, estimated within frogs; population differences in load and detection are absorbed by the frog-level intercepts. With `collect_lab` included (and pooled `sigma` and swab SD), estimates were: PrepMan −0.38 log10 (95% CI −0.65 to −0.10); SNARL qPCR +0.06 log10 (−0.08 to +0.20); `s` 0.05 (−0.04 to 0.14); qPCR noise 0.28 and swab noise 0.40 log10; no divergences.
   - **Model 2b (constant bias)**: Model 2 without `s`. Compared with Model 2 by PSIS-LOO; if comparable (`elpd_diff` < ~2 SE), the constant-bias model is preferred and the correction to the Qiagen–RZ scale is a single constant per pipeline (exported to `fits/m2b_constant_correction_export.rds`). Expected from earlier fits: SNARL/PrepMan → RZ/Qiagen ≈ +0.32 log10 (≈ 2.1×).
   - **LOO result**: `elpd_diff` = −0.6 (SE 2.5) for Model 2b vs Model 2, i.e. no detectable difference in predictive performance. The slope adds nothing, so **Model 2b (constant bias) is preferred**: there is no evidence of a load-dependent SNARL qPCR bias.
   - **No extraction × qPCR interaction**: all four pipelines are observed, but each pair compares different frog populations (see panel table), so an interaction would be confounded with population.
   - **First-version results** (lab-specific `sigma`): PrepMan reads ≈ 0.50 log10 lower than Qiagen (95% CI 0.20–0.78) and has more non-detects; no evidence of a SNARL qPCR bias (≈ +0.10 log10, CI −0.05 to +0.26) or slope (`s` = 0.05, CI −0.04 to 0.14); SNARL qPCR has borderline more non-detects (`hu` +1.06, CI −0.04 to 2.20). The hypothesis is partly supported: SNARL qPCR does not under-quantify positives but may miss more low-load samples.
   - **Correction to a common scale**: positive readings from any pipeline are converted to the reference (Qiagen + RZ) scale by inverting the model, `L = (y − d_ext − d_q + s·L_ref) / (1 + s)`; correction parameter draws and `L_ref` are exported to `fits/m2_correction_export.rds`. Zeros cannot be corrected.
   - **Priors from Model 1** (SDs widened ×2) only where Model 1 is well identified and parameters mean the same thing: `sigma` (qPCR replicate noise, pooled), swab SD in `L` (√(swab noise² − qPCR noise²)), frog-level `hu` SD. Intercepts get weakly informative priors. The `L`–`hu` correlation gets `lkj(2)`, since brms does not support centred priors on correlations.

### Key questions
1) Are between-replicate quantitative differences significantly different between labs?
2) Are quantitative outputs from the qPCR step significantly different between labs? If not, how do they differ, in the interest of making the results comparable?
3) Are quantitative outputs from the extraction step significantly different between labs? If not, how do they differ, in the interest of making the results comparable?

## Files
Rendered output (`*.html`, `*_files/`) and fitted models are not tracked in git. Fitted models are saved to `$data_dir/bd_qpcr_results/rz_snarl_qpcr_comparison/fits/`.

- `R/data_import.R`: pulls experiment data (xlsx) and database results, and saves them to `rz_snarl_qpcr_raw_data.RData` in the data directory. Run once, and again whenever source data change (only step needing a database connection).
- `R/data_prep.R`: loads the saved raw data and builds cleaned experiment (`snarl_clean`, `snarl_wide`) and replicate (`replicates_clean`, etc.) tables. Sourced by all qmd files.
- `rz_snarl_qpcr_diagnostics.qmd`: data exploration and diagnostics.
  - **Experiment data**: Bd load by frog, wide across extraction method × qPCR lab.
  - **Data Summary**: pairwise scatter of experiment Bd loads across extraction method × qPCR lab (panels A–E).
  - **Replicate Variability**: all within-group replicate pairs, faceted by qPCR lab × replicate type, kept separate from the experiment data. Used to assess whether replicate variability can be pooled across labs.
  - **Replicate Diagnostics**: pair agreement (detection outcomes, SD of log10 differences) for experiment panels vs replicates; swab replicate group checks (sample types, name conflicts, plates, projects); largest swab disagreements; difference vs. mean (Bland–Altman) plots to test whether noise depends on load and population.
- `rz_snarl_qpcr_model1_replicates.qmd`: Model 1 (replicates). Model data, structure (family, formula, priors), fit, convergence & correlations, detection & noise, posterior predictive checks. Exports posterior draws and summaries of population-level effects, group-level SDs, and correlations to `fits/m1_replicates_export.rds`.
- `rz_snarl_qpcr_model2_experiment.qmd`: Model 2 (experiment). Reads the Model 1 export (render Model 1 first; warns if Model 1 was fit on an older data import). Model data, priors from Model 1, structure, fit, summary, lab differences table (log10 biases, slope, detection, noise), divergence pairs plot, SNARL − RZ qPCR difference across load and detection by pipeline, posterior predictive checks, Model 2b (constant bias) fit and LOO comparison with constant correction table, correction to the Qiagen–RZ scale under Model 2 (plot, table, exported draws).

In both model files, set `refit = TRUE` in the setup chunk to force a refit; otherwise saved fits are loaded (brms refits automatically if the formula or data change).
