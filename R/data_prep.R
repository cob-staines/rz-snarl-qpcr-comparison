# Data prep for RZ-SNARL qPCR comparison
# Loads raw data saved by R/data_import.R (no database connection needed), and builds:
#   snarl_clean, snarl_wide         experiment data, long & wide by frog
#   db_clean                        database ITS1 results (experiment swabs excluded)
#   qpcr_replicates_clean, swab_replicates_clean, replicates_clean, replicates_wide
# Sourced by rz_snarl_qpcr_diagnostics.qmd and rz_snarl_qpcr_modeling.qmd

librarian::shelf(tidyverse, here, janitor)

wddir = here(Sys.getenv("data_dir"), "bd_qpcr_results", "rz_snarl_qpcr_comparison")

# loads snarl_raw, bd_results, import_date
load(here(wddir, "rz_snarl_qpcr_raw_data.RData"))
message("Using data imported on ", import_date)

snarl_clean = snarl_raw %>%
  clean_names() %>%
  rename(extraction_lab = lab_extract,
         extraction_method = extract_method,
         qpcr_lab = lab_qpcr,
         bd_quant_cycle = quant_cycle,
         bd_start_quant = start_quant,
         frog_id = frog_num) %>%
  mutate(extraction_lab = recode_values(extraction_lab,
                                        "pitt" ~ "rz",
                                        "snarl" ~ "snarl",
                                        default = NA),
         collect_lab = recode_values(collect_lab,
                                     "pitt" ~ "rz",
                                     "snarl" ~ "snarl",
                                     default = NA),
         qpcr_lab = recode_values(qpcr_lab,
                                  "pitt" ~ "rz",
                                  "snarl" ~ "snarl",
                                  default = NA),
         frog_id = as.character(frog_id),
         src = "snarl")

# catch unexpected lab values silently recoded to NA
stopifnot(!anyNA(snarl_clean$extraction_lab),
          !anyNA(snarl_clean$collect_lab),
          !anyNA(snarl_clean$qpcr_lab))

snarl_wide = snarl_clean %>%
  pivot_wider(id_cols = c("frog_id", "collect_lab"),
              names_from = c("extraction_method", "qpcr_lab"),
              values_from = "bd_load")

# database results, excluding experiment swabs so they are not double-counted as replicates
db_clean = bd_results %>%
  filter(!sample_name_bd %in% snarl_clean$bd_swab_id) %>%
  rename(bd_swab_id = sample_name_bd,
         bd_quant_cycle = bd_cycle_quant,
         bd_start_quant = bd_target_quant,
         bd_load = bd_its1_copies_per_swab,
         notes = comments_qpcr) %>%
  mutate(extraction_lab = recode_values(extraction_lab,
                                        "rz_lab" ~ "rz",
                                        default = extraction_lab),
         qpcr_lab = recode_values(qpcr_lab,
                                  "rz_lab" ~ "rz",
                                  default = qpcr_lab),
         extraction_method = recode_values(extraction_kit,
                                           "qiagen_dneasy" ~ "qiagen",
                                           default = extraction_kit),
         collect_lab = NA_character_,  # not recorded in db; not used for replicates
         frog_id = capture_id,
         src = "db",
         species_capture = str_to_sentence(gsub("_", " ", taxon_capture)))

# qPCR replicates: same extract (swab + extraction plate) run more than once within a qPCR lab
qpcr_replicates_clean = db_clean %>%
  group_by(bd_swab_id, extraction_plate_name, qpcr_lab) %>%
  filter(n() > 1) %>%
  mutate(replicate_type = "qpcr",
         replicate_group = paste(bd_swab_id, extraction_plate_name, qpcr_lab, sep = "_"),
         replicate_id = row_number()) %>%
  ungroup()

# swab replicates: multiple swabs from the same capture, through the same extraction & qPCR lab
swab_replicates_clean = db_clean %>%
  # keep one result per swab (first qPCR plate), so qPCR re-runs are not counted as swab replicates
  arrange(qpcr_plate_name) %>%
  group_by(bd_swab_id, extraction_lab, extraction_method, qpcr_lab) %>%
  slice_head(n = 1) %>%
  group_by(capture_id, extraction_lab, extraction_method, qpcr_lab) %>%
  filter(n() > 1) %>%
  mutate(replicate_type = "swab",
         replicate_group = paste(capture_id, extraction_lab, qpcr_lab, sep = "_"),
         replicate_id = row_number()) %>%
  ungroup()

replicates_clean = bind_rows(qpcr_replicates_clean,
                             swab_replicates_clean) %>%
  select(all_of(colnames(snarl_clean)),
         qpcr_plate_name,
         project_id,
         replicate_type,
         replicate_group,
         replicate_id)

replicates_wide = replicates_clean %>%
  pivot_wider(id_cols = c("replicate_type", "replicate_group", "frog_id",
                          "extraction_method", "qpcr_lab"),
              names_from = "replicate_id",
              values_from = "bd_load")
