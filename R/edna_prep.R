# eDNA data prep for RZ-SNARL qPCR comparison
# Loads RZ eDNA plate data saved by R/edna_import.R and builds:
#   edna_clean           one row per usable Bd well (latest extraction & run per sample)
#   edna_units           one row per sample (its wells on its latest run)
#   edna_cleaning_log    wells & samples remaining after each cleaning step
#   edna_excluded_plates plates dropped for a positive negative control
# Sourced by the diagnostics and Model 0 qmd files

librarian::shelf(tidyverse, here)

wddir = here(Sys.getenv("data_dir"), "bd_qpcr_results", "rz_snarl_qpcr_comparison")

# loads edna_wells, edna_plates, edna_import_date
load(here(wddir, "rz_edna_raw_data.RData"))
message("Using eDNA data imported on ", edna_import_date)

# eDNA quantities are copies per reaction; x40 puts them on the experiment's bd_load scale
# (ITS1 copies per swab, Qiagen extracts: bd_load = 40 x bd_start_quant)
edna_load_factor = 40

# Ct range of ~28 between wells (one well amplified at Ct < 10)
edna_excluded_samples = c("PE_24_EDNA_205")

edna_named = edna_wells %>%
  # unnamed wells are empty, with spurious amplification at very low Ct
  filter(!is.na(sample_name)) %>%
  mutate(name_std = sample_name %>%
           str_squish() %>%
           str_to_upper() %>%
           str_replace_all("[ -]+", "_") %>%
           str_replace_all("(?<=_)0+(?=\\d)", ""),  # PAN_22_054 = PAN_22_54
         is_ntc = task == "NTC",
         is_field_nc = str_detect(name_std, "(^|_)NC(_?\\d+)?(_|$)"),
         detected = !is.na(bd_ct))

# plates with any positive no-template or field negative control
edna_excluded_plates = edna_named %>%
  filter(is_ntc | is_field_nc) %>%
  group_by(plate, region) %>%
  summarise(positive_controls = paste(sample_name[detected], collapse = ", "),
            .groups = "drop") %>%
  filter(positive_controls != "")

edna_unknowns = edna_named %>%
  filter(task == "UNKNOWN", !is_field_nc) %>%
  # 2022-era plates label triplicate wells _A/_B/_C (one well each, adjacent, same plate)
  mutate(well_letter = str_extract(name_std, "(?<=_)[ABC]$"),
         name_stem = if_else(is.na(well_letter), name_std, str_remove(name_std, "_[ABC]$"))) %>%
  group_by(plate, name_stem) %>%
  mutate(letters_are_wells = !anyNA(well_letter) & n() >= 2 & n_distinct(well_letter) == n()) %>%
  ungroup() %>%
  mutate(name_extract = if_else(letters_are_wells, name_stem, name_std),
         # remaining B/C suffixes (e.g. Pe_25_511_B, Pe 24 258B) are re-extractions of the same sample
         extraction = coalesce(str_remove(str_extract(name_extract, "(?<=\\d)_?[BC]$"), "_"), "A"),
         sample = str_remove(name_extract, "(?<=\\d)_?[BC]$"),
         is_tissue = str_detect(sample, "^TISSUE"))

log_step = function(df, step) {
  tibble(step = step,
         wells = nrow(df),
         samples = n_distinct(df$sample),
         samples_detected = n_distinct(df$sample[df$detected]))
}

# re-runs and re-extractions were done when contamination was suspected:
# keep the latest extraction of each sample, then its latest run
edna_latest = edna_unknowns %>%
  group_by(sample) %>%
  filter(extraction == max(extraction)) %>%
  filter(run_time == max(run_time)) %>%
  ungroup()

edna_no_contam = edna_latest %>%
  filter(!plate %in% edna_excluded_plates$plate)

# IPC not amplified = inhibition (matches the lab's Working_Results IPC calls)
edna_ipc = edna_no_contam %>%
  filter(!is.na(ipc_ct))

edna_clean = edna_ipc %>%
  filter(!sample %in% edna_excluded_samples) %>%
  transmute(sample, region, plate, run_time, well, sample_name, extraction, is_tissue,
            bd_ct,
            bd_start_quant = if_else(detected, bd_quantity, 0),
            bd_load = edna_load_factor * bd_start_quant,
            detected)

stopifnot(!anyNA(edna_clean$bd_load))

edna_cleaning_log = bind_rows(
  log_step(edna_unknowns, "named sample wells (no controls)"),
  log_step(edna_latest, "latest extraction & run per sample"),
  log_step(edna_no_contam, "drop plates with positive negative control"),
  log_step(edna_ipc, "drop IPC-failed wells"),
  log_step(edna_clean, "drop excluded samples"))

edna_units = edna_clean %>%
  group_by(sample, region, plate, is_tissue) %>%
  summarise(n_wells = n(),
            n_detected = sum(detected),
            mean_log10_load = if_else(n_detected > 0, mean(log10(bd_load[detected])), NA_real_),
            sd_log10_load = if_else(n_detected > 1, sd(log10(bd_load[detected])), NA_real_),
            .groups = "drop")
