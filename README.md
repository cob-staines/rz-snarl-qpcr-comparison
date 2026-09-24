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
Bayesian hurdle-lognormal models in brms (cmdstanr backend) on Bd load (ITS1 copies per swab). Each unit (replicate group, or frog) has a latent load expressed as a group-level intercept that is correlated across `mu`, `hu`, and `sigma`; the correlations capture load-dependent detection (zeros treated as non-detections) and load-dependent noise.

Bd load within a population is generally modeled as a hurdle-lognormal quantity. Here, the goal is not to represent a population, but to compare relative values between control groups. The model must therefore handle zeros as non-detections (possibly load-dependent detection), rather than as true absence.

Two-stage approach:
1. **Model 1: Replicates**. Correlated group-level intercepts per replicate group in `mu`, `hu`, `sigma`; qPCR lab effect on `hu`; `sigma` by noise type (swab replicates: total noise; qPCR replicates: within-plate well noise); qPCR run (extract × plate) effect in `mu`. Posterior summaries are saved as priors for Model 2.
2. **Model 2: Experiment** (to be implemented). Latent load per frog; lab-specific extraction and qPCR effects on `mu` and `hu`; priors on noise and detection informed by Model 1.

### Key questions
1) Are between-replicate quantitative differences significantly different between labs?
2) Are quantitative outputs from the qPCR step significantly different between labs? If not, how do they differ, in the interest of making the results comparable?
3) Are quantitative outputs from the extraction step significantly different between labs? If not, how do they differ, in the interest of making the results comparable?

## Files
Rendered output (`*.html`, `*_files/`) and fitted models are not tracked in git. Fitted models are saved to `$data_dir/bd_qpcr_results/rz_snarl_qpcr_comparison/fits/`.

- `R/data_prep.R`: loads experiment data (xlsx) and database results, and builds cleaned experiment (`snarl_clean`, `snarl_wide`) and replicate (`replicates_clean`, etc.) tables. Sourced by both qmd files.
- `rz_snarl_qpcr_diagnostics.qmd`: data exploration and diagnostics.
  - **Experiment data**: Bd load by frog, wide across extraction method × qPCR lab.
  - **Data Summary**: pairwise scatter of experiment Bd loads across extraction method × qPCR lab (panels A–E).
  - **Replicate Variability**: all within-group replicate pairs, faceted by qPCR lab × replicate type, kept separate from the experiment data. Used to assess whether replicate variability can be pooled across labs.
  - **Replicate Diagnostics**: pair agreement (detection outcomes, SD of log10 differences) for experiment panels vs replicates; swab replicate group checks (sample types, name conflicts, plates, projects); largest swab disagreements; difference vs. mean (Bland–Altman) plots to test whether noise depends on load and population.
- `rz_snarl_qpcr_modeling.qmd`: model fitting and results. Set `refit = TRUE` in the setup chunk to refit; otherwise saved fits are loaded.
  - **Model 1: Replicates**: model data, structure (family, formula, priors), fit, convergence & correlations, detection & noise vs. latent load, posterior predictive checks, prior summaries for Model 2.
  - **Model 2: Experiment**: placeholder.
