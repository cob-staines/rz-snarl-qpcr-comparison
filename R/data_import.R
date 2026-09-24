# Data import for RZ-SNARL qPCR comparison
# Pulls raw experiment data (xlsx) and RIBBiTR database results, and saves them to
#   $data_dir/bd_qpcr_results/rz_snarl_qpcr_comparison/rz_snarl_qpcr_raw_data.RData
# Run once (and whenever the source data change); R/data_prep.R loads the saved file.

librarian::shelf(tidyverse, dbplyr, here, rio, RPostgres, DBI, RIBBiTR-BII/ribbitrrr)

wddir = here(Sys.getenv("data_dir"), "bd_qpcr_results", "rz_snarl_qpcr_comparison")

# establish database connection
dbcon = hopToDB("railway", hopReg = FALSE)

# load data
snarl_raw = import(here(wddir, "qpcr_compare_data_260904.xlsx"),
                   na = c("NA", "NaN"))

# results from database
db_bd = tbl(dbcon, Id("survey_data", "bd_qpcr_results"))
db_sample = tbl(dbcon, Id("survey_data", "sample"))
db_capture = tbl(dbcon, Id("survey_data", "capture"))
db_survey = tbl(dbcon, Id("survey_data", "survey"))
db_visit = tbl(dbcon, Id("survey_data", "visit"))
db_site = tbl(dbcon, Id("survey_data", "site"))

# filter before counting replicates, so replicate counts reflect usable results only
bd_results = db_bd %>%
  inner_join(db_sample, by = "sample_id") %>%
  inner_join(db_capture, by = "capture_id") %>%
  left_join(db_survey %>% select(survey_id, visit_id), by = "survey_id") %>%
  left_join(db_visit %>% select(visit_id, date, site_id, project_id), by = "visit_id") %>%
  left_join(db_site %>% select(site_id, site), by = "site_id") %>%
  filter(ipc_pass | is.na(ipc_pass),
         standard_target_type == "ITS1",
         !is.na(bd_its1_copies_per_swab),
         qpcr_lab %in% c("snarl", "rz_lab")) %>%
  select(result_id,
         sample_id,
         capture_id,
         sample_name_bd,
         bd_cycle_quant,
         bd_target_quant,
         bd_its1_copies_per_swab,
         extraction_plate_name,
         extraction_lab,
         qpcr_plate_name,
         qpcr_lab,
         extraction_kit,
         taxon_capture,
         comments_qpcr,
         # metadata for replicate diagnostics
         sample_type,
         sample_name_conflict,
         swab_type,
         qpcr_well_replicates = replicates,
         extraction_date,
         qpcr_date,
         visit_date = date,
         site,
         project_id) %>%
  collect()

dbDisconnect(dbcon)

import_date = Sys.Date()

save(snarl_raw, bd_results, import_date,
     file = here(wddir, "rz_snarl_qpcr_raw_data.RData"))
