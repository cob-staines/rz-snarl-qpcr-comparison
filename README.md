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

## Hypothesis
qPCR results from the SNARL lab underestimate low Bd quantities compared to the RZ qPCR protocol, but are otherwise comparable.

## Methods
Using a Bayesian GAM model (brms) to model the expected transformation function of the extraction and qPCR steps at each lab, taking into consideration the noise term (lognormal).

Bd load within a population is generally modeled as a hurdle-lognormal quantity. Here, the goal is not to represent a population, but to compare relative values between control groups. The model must therefore handle zeros as non-detections (possibly load-dependent detection), rather than as true absence.

### Key questions
1) Are between-replicate quantitative differences significantly different between labs?
2) Are quantitative outputs from the qPCR step significantly different between labs? If not, how do they differ, in the interest of making the results comparable?
3) Are quantitative outputs from the extraction step significantly different between labs? If not, how do they differ, in the interest of making the results comparable?

## Files
- `rz_snarl_qpcr_comparison.qmd`: currently the key data wrangling, modeling, and reporting file. We might split this into separate files if it becomes too clunky. Rendered output (`*.html`, `*_files/`) is gitignored.
