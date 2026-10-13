# eDNA data import for RZ-SNARL qPCR comparison
# Reads RZ lab eDNA Bd qPCR plate exports (triplicate wells of one extract) and saves them to
#   $data_dir/bd_qpcr_results/rz_snarl_qpcr_comparison/rz_edna_raw_data.RData
# Run once (and whenever new plates are added); R/edna_prep.R loads the saved file.

librarian::shelf(tidyverse, here, readxl)

wddir = here(Sys.getenv("data_dir"), "bd_qpcr_results", "rz_snarl_qpcr_comparison")
edna_dir = here(Sys.getenv("data_dir"), "edna_results", "edna_bd_qpcr", "Complete")

edna_files = list.files(edna_dir, pattern = "\\.xlsx$", recursive = TRUE, full.names = TRUE)

edna_files = tibble(path = edna_files,
                    plate = basename(edna_files),
                    region = basename(dirname(edna_files))) %>%
  # macOS resource-fork files
  filter(!str_detect(plate, "^\\._"),
         # frog swabs, not eDNA
         !str_detect(plate, "^Lehtinen_Swabs")) %>%
  # some plates are filed under more than one region (e.g. PennPAN2022 plate3and1)
  arrange(plate, region) %>%
  distinct(plate, .keep_all = TRUE)

# one row per well & target from the "Results" sheet; the header row position varies between exports
read_edna_plate = function(path, plate, region) {
  raw = read_excel(path, sheet = "Results", col_names = FALSE, col_types = "text", .name_repair = "minimal")
  header_row = which(raw[[1]] == "Well")[1]
  run_end = raw[[2]][which(raw[[1]] == "Experiment Run End Time")[1]]

  read_excel(path, sheet = "Results", skip = header_row - 1, col_types = "text", .name_repair = "minimal") %>%
    transmute(plate = plate,
              region = region,
              run_end = run_end,
              well = `Well Position`,
              sample_name = `Sample Name`,
              target = `Target Name`,
              task = Task,
              ct = suppressWarnings(as.numeric(CT)),             # "Undetermined" -> NA
              quantity = suppressWarnings(as.numeric(Quantity)), # copies per reaction
              curve_slope = suppressWarnings(as.numeric(Slope)),
              curve_efficiency = suppressWarnings(as.numeric(Efficiency)),
              curve_r2 = suppressWarnings(as.numeric(`R(superscript 2)`)))
}

edna_results_long = pmap_dfr(edna_files %>% select(path, plate, region), read_edna_plate) %>%
  filter(!is.na(well))

# wide: one row per well, Bd quantity and IPC Ct side by side
edna_wells = edna_results_long %>%
  filter(target %in% c("Bd", "IPC")) %>%
  select(plate, region, run_end, well, sample_name, target, task, ct, quantity) %>%
  pivot_wider(names_from = target, values_from = c(task, ct, quantity)) %>%
  transmute(plate, region, well, sample_name,
            task = task_Bd,
            bd_ct = ct_Bd,
            bd_quantity = quantity_Bd,
            ipc_ct = ct_IPC,
            run_time = parse_date_time(str_remove(run_end, " [A-Z]{3,4}$"), "Ymd HMS p", quiet = TRUE))

# Bd standard curve per plate
edna_plates = edna_results_long %>%
  filter(target == "Bd") %>%
  group_by(plate, region) %>%
  summarise(curve_slope = first(na.omit(curve_slope)),
            curve_efficiency = first(na.omit(curve_efficiency)),
            curve_r2 = first(na.omit(curve_r2)),
            .groups = "drop") %>%
  left_join(edna_wells %>% distinct(plate, run_time), by = "plate")

stopifnot(!anyNA(edna_plates$run_time))

edna_import_date = Sys.Date()

save(edna_wells, edna_plates, edna_import_date,
     file = here(wddir, "rz_edna_raw_data.RData"))
