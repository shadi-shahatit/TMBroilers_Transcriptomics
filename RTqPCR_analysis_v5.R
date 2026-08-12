#### RT-qPCR analysis - v5
#### Shadi Shahatit, MBV AlZghoul Lab, JUST, 2025

# This script reads raw qPCR Excel sheets, containing sample ID, housekeeping genes, and multiple target-gene Ct values.
# It standardizes sample names, extracts biological replicates, reshapes the data 
# from wide to long format, assigns technical replicates when present, and computes 
# ΔCt, ΔΔCt, and fold-change values for every gene under every condition using both 
# ACTIN and GAPDH. For each sheet and each housekeeping gene, it generates replicate-level
# data, biological summaries, significance testing (t-test or Wilcoxon), and fold-change plots, 
# and finally exports clean per-gene tables (Ct1/Ct2, ΔCt, ΔΔCt, fold-change, and sum stats) in Excel.

# INPUT:    Excel file of Ct values (sheets = primer concentration; rows = samples; columns = genes).
# OUTPUT:   Dataframes of Ct tables with replicate and biological-level data, including ΔCt, ΔΔCt, fold-change summaries.
# OUTPUT:   Bar plots of fold change in reference to the control.
# OUTPUT:   Excel file containing per-gene tables (Ct1/Ct2, ΔCt, ΔΔCt, fold-change, and summary statistics).

# USER-DEFINED PARAMETERS: Excel file path; Excel column headers; condition map; calibrator condition label 

## Load the packages  ======================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(rstatix)
library(tibble)
library(writexl)

## Define your parameters  ======================

## path for your CT raw data form the machine
## wide-sheet settings (sheets are primer concentrations; rows are samples; columns are genes)

#### TM_chicken 

## for MultTissue_Devo
# qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_MultTissue_Devo\\Liver\\OPT_LIVER_DEVO_SS_100_mod.xlsx"
# qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_MultTissue_Devo\\Muscle\\OPT_MUSCLE_DEVO_SS_100_mod.xlsx"
# qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_MultTissue_Devo\\Spleen\\OPT_SPLEEN_DEVO_SS_100_mod.xlsx"

qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_MultTissue_Devo\\Liver\\Liver_Devo_CT1_CT2_cDNA_Final_SS_mod.xlsx"
# qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_MultTissue_Devo\\Muscle\\Muscle_Devo_CT1_CT2_cDNA_Final_SS_mod.xlsx"
# qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_MultTissue_Devo\\Spleen\\Spleen_Devo_CT1_CT2_cDNA_Final_SS_mod.xlsx"

## for muscle
# qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_muscle\\MUSCLE_CT1_CT2_cDNA_1_400_5_2_5_mod_SS.xlsx"
## for liver
# qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_liver\\LIVER_CT1_CT2_cDNA_1_XXX_5_2_5_SS.xlsx"
## for spleen
# qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_spleen\\CT1_CT2_cDNA_1_200_mod_SS.xlsx"
# qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_spleen\\CT1_CT2_cDNA_1_400_mod_SS.xlsx"

#### Atro_Ex

# ## for Atro_Ex
# qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_Atrophy\\CT1_CT2_cDNA_1_400_primers_5x_mod_SS.xlsx"
# qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_Atrophy\\OPT_1_100_cDNA_10X_5X_2.5X_PRIMERS_mod.xlsx"
# qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_Atrophy\\OPT_1_200_cDNA_10X_5X_2.5X_PRIMERS_mod.xlsx"
# qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_Atrophy\\OPT_1_400_cDNA_10X_5X_2.5X_PRIMERS_mod.xlsx"

# qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_Atrophy\\Heart\\Atrophy_Ex_project_heart_CT1_CT2_24_SAMPLES_FINAL_SS_mod.xlsx"

# qpcr_file_path <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_Atrophy\\Muscle\\Atrophy_Ex_project_muscle_CT1_CT2_FINAL_SS_mod.xlsx"

## list the sheets to read; NULL = read all (e.g., "concentration a", "concentration b")
sheets_to_use  <- NULL

## exact headers
sample_id_col  <- "SAMPLE ID"
actin_col      <- "B-ACTIN"
gapdh_col      <- "GAPDH"

## remove "5X" or "2.5X" tokens from SAMPLE ID; the concentration already exists in the sheet name
strip_leading_token <- TRUE

## targets = all non-ID/HK by default; set a character vector to force a subset
target_gene_cols <- NULL

## condition map

## Order matters: put longer/more specific patterns first (e.g. HLEM before HLM) for grepl

# condition_map <- c(
#   "D19L1"   = "D19Liver",
#   "D7L2"    = "D7Liver",
#   "D22L2" = "D22Liver"
#   )

# condition_map <- c(
#   "D19M2"   = "D19Muscle",
#   "D7M1"    = "D7Muscle",
#   "D22M1" = "D22Muscle"
# )

# condition_map <- c(
#   "D19S1"   = "D19Spleen",
#   "D7S2"    = "D7Spleen",
#   "D22S1" = "D22Spleen"
# )

condition_map <- c(
  ## Liver
  "D19L1"   = "D19Liver",
  "D19L2"   = "D19Liver",

  "D7L2"    = "D7Liver",
  "D7L3"    = "D7Liver",

  "D22L2" = "D22Liver",
  "D22L3" = "D22Liver",
  "D22L4" = "D22Liver"
    
  # ## Muscle
  # "D19M2"   = "D19Muscle",
  # "D19M3"   = "D19Muscle",
  # 
  # "D7M1"    = "D7Muscle",
  # "D7M3"    = "D7Muscle",
  # 
  # "D22M1" = "D22Muscle",
  # "D22M2" = "D22Muscle",
  # "D22M3" = "D22Muscle"

  # ## Spleen
  # "D19S1"   = "D19Spleen",
  # "D19S3"   = "D19Spleen",
  # 
  # "D7S2"    = "D7Spleen",
  # "D7S3"    = "D7Spleen",
  # 
  # "D22S1" = "D22Spleen",
  # "D22S2" = "D22Spleen",
  # "D22S6" = "D22Spleen"
)

# condition_map <- c(
#   "CON_N"   = "Con_N",
#   "TM_N"    = "TM_N",
#   "CON_AHS" = "Con_AHS",
#   "TM_AHS"  = "TM_AHS",
#   "CON_CHS" = "Con_CHS",
#   "TM_CHS"  = "TM_CHS"
#   )

## control condition - could change if other comp is needed

calibrator_condition_label <- "D19Liver"
# calibrator_condition_label <- "D19Muscle"
# calibrator_condition_label <- "D19Spleen"

# calibrator_condition_label <- "Con_N"

## desired condition levels

## for MultTissue_Devo
condition_levels <- c("D19Liver","D7Liver","D22Liver")
# condition_levels <- c("D19Muscle","D7Muscle","D22Muscle")
# condition_levels <- c("D19Spleen","D22Spleen","D7Spleen")

## for muscle
# condition_levels <- c("Con_N","TM_N","Con_AHS","TM_AHS","Con_CHS","TM_CHS")
## for liver
# condition_levels <- c("Con_N","TM_N","Con_AHS","TM_AHS")
## for spleen
# condition_levels <- c("Con_N","TM_N","Con_AHS","TM_AHS","Con_CHS","TM_CHS")

#### Atro_Ex - Muscle

# condition_map <- c(
#   "HLEM" = "HLEM",
#   "HLM"  = "HLM",
#   "HEM"  = "HEM",
#   "LEM"  = "LEM",
#   "CM"   = "CM"
#   )

## control condition - could change if other comp is needed

# calibrator_condition_label <- "CM"

## desired condition levels

## for Atro_Ex - Muscle
# condition_levels <- c("CM","HLM","HLEM","HEM","LEM")

## for Atro_Ex - Intestine

# condition_map <- c(
#   "HLET" = "HLET",
#   "HLT"  = "HLT",
#   "HET"  = "HET",
#   "LET"  = "LET",
#   "CT"   = "CT"
# )

## control condition - could change if other comp is needed

# calibrator_condition_label <- "CT"

## desired condition levels

## for Atro_Ex - Intestine
# condition_levels <- c("CT","HLT","HLET","HET","LET")

## for Atro_Ex - Heart

# condition_map <- c(
#   "HLEH" = "HLEH",
#   "HLH"  = "HLH",
#   "HEH"  = "HEH",
#   "LEH"  = "LEH",
#   "CH"   = "CH"
# )

## control condition - could change if other comp is needed

# calibrator_condition_label <- "CH"

## desired condition levels

## for Atro_Ex - Heart
# condition_levels <- c("CH","HLH","HLEH","HEH","LEH")

## Read target gene sheets  ======================

all_sheet_names <- readxl::excel_sheets(qpcr_file_path)
if (is.null(sheets_to_use)) sheets_to_use <- all_sheet_names
if (!length(sheets_to_use)) stop("No sheets to process in the file.")

## helpers

## read_one: a function to read one Excel sheet and tag rows with the sheet name
read_one <- function(sh) {
  df <- readxl::read_excel(qpcr_file_path, sheet = sh)
  df$.__sheet__ <- sh
  df
}

## sanitize_sheet_name: a function to convert a sheet name into a safe identifier
sanitize_sheet_name <- function(x) {
  x %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("[^A-Za-z0-9_]+", "_")
}

## Validate and prepare per-sheet data  ======================

raw_list <- purrr::map(sheets_to_use, read_one)
names(raw_list) <- sheets_to_use

## ensure required columns exist per sheet, determine target genes, and force numeric
raw_list <- purrr::imap(raw_list, function(df, sh) {
  if (!(sample_id_col %in% names(df))) {
    stop("Required SAMPLE ID column not found in sheet '", sh, "': ", sample_id_col)
  }
  if (!all(c(actin_col, gapdh_col) %in% names(df))) {
    stop("Both housekeeping columns must exist in sheet '", sh, "': ", actin_col, " and ", gapdh_col)
  }
  default_exclude <- c(sample_id_col, actin_col, gapdh_col, ".__sheet__")
  tg <- if (is.null(target_gene_cols)) setdiff(names(df), default_exclude) else target_gene_cols
  if (!length(tg)) stop("No target gene columns detected in sheet '", sh, "'.")
  num_cols <- intersect(c(actin_col, gapdh_col, tg), names(df))
  df <- df %>% mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.))))
  attr(df, "target_gene_cols") <- tg
  df
})

## Build replicate-level target CT table for all sheets  ======================

## parse_sample_id: a function to strip concentration tokens, normalize sample IDs, match conditions via condition_map, and extract biological replicate numbers
parse_sample_id <- function(sample_ids, strip_leading_token = TRUE) {
  
  purrr::map_dfr(sample_ids, function(x) {
    x <- as.character(x) %>% stringr::str_squish()
    
    # strip leading concentration token e.g. "10X ", "2.5X "
    if (strip_leading_token) {
      x <- sub("^[0-9]+(?:\\.[0-9]+)?\\s*[xX]\\s*", "", x, perl = TRUE)
    }
    
    # normalize to uppercase, collapse separators
    x_norm <- x %>%
      toupper() %>%
      gsub("[-\\s]+", "_", x = ., perl = TRUE) %>%
      stringr::str_squish()
    
    # match condition via grepl on the map
    matched_label <- NA_character_
    for (pattern in names(condition_map)) {
      if (grepl(pattern, x_norm, ignore.case = TRUE)) {
        matched_label <- condition_map[[pattern]]
        break
      }
    }
    
    # extract trailing replicate number if present
    rep_num <- stringr::str_extract(x_norm, "\\d+$") %>% # replacing str_match(sample_ids, "-(\\d+)$")
      suppressWarnings(as.integer())
    tibble::tibble(condition_label = matched_label, bio_replicate = rep_num)
  }) %>%
    dplyr::group_by(condition_label) %>%
    dplyr::mutate(
      bio_replicate = ifelse(is.na(bio_replicate), row_number(), bio_replicate)
    ) %>%
    dplyr::ungroup()
}

## build_base_long_unified: a function to attach parsed sample IDs, attach housekeeping Cts, pivot wide to long, and apply condition ordering
build_base_long_unified <- function(
    df,
    sample_id_col,
    actin_col,
    gapdh_col,
    target_gene_cols  = NULL,
    strip_leading_token = TRUE
) {
  
  # infer target genes if not provided
  default_exclude <- c(sample_id_col, actin_col, gapdh_col, ".__sheet__")
  tg <- if (is.null(target_gene_cols)) setdiff(names(df), default_exclude) else target_gene_cols
  
  # parse sample IDs into condition_label + bio_replicate
  id_parsed <- parse_sample_id(df[[sample_id_col]], strip_leading_token = strip_leading_token)
  
  # attach parsed IDs and housekeeping Cts
  df2 <- bind_cols(df, id_parsed) %>%
    group_by(condition_label) %>%
    # mutate(bio_replicate = ifelse(is.na(bio_replicate), row_number(), bio_replicate)) %>%
    ungroup() %>%
    mutate(
      ct_actin = .[[actin_col]],
      ct_gapdh = .[[gapdh_col]]
    )
  
  # pivot to long, drop missing Ct_target
  long <- df2 %>%
    select(.__sheet__, all_of(sample_id_col), condition_label, bio_replicate,
           ct_actin, ct_gapdh, all_of(tg)) %>%
    pivot_longer(
      cols      = all_of(tg),
      names_to  = "target_gene",
      values_to = "ct_target"
    ) %>%
    filter(!is.na(ct_target)) %>%
    mutate(condition_label = as.character(condition_label) %>% str_squish())
  
  # apply condition ordering from global condition_levels
  present <- unique(long$condition_label)
  ordered <- c(condition_levels[condition_levels %in% present], setdiff(present, condition_levels))
  long$condition_label <- factor(long$condition_label, levels = ordered)
  
  # assign technical replicate index per (gene, condition, bio) if duplicates exist
  long <- long %>%
    group_by(target_gene, condition_label, bio_replicate) %>%
    arrange(.data[[sample_id_col]], .by_group = TRUE) %>%
    mutate(tech_replicate = row_number()) %>%
    ungroup()
  
  long
}

## 2 remove

# ## build_base_long_unified: a function to parse SAMPLE ID, standardize conditions, attach housekeeping Ct, pivot wide→long, and apply condition ordering
# build_base_long_unified <- function(
#     df,
#     sample_id_col,
#     actin_col,
#     gapdh_col,
#     target_gene_cols = NULL,
#     # for muscle
#     # condition_levels = c("Con_N","TM_N","Con_AHS","TM_AHS","Con_CHS","TM_CHS"),
#     # for liver
#     # condition_levels = c("Con_N","TM_N","Con_AHS","TM_AHS"),
#     # for spleen
#     # condition_levels = c("Con_N","TM_N","Con_AHS","TM_AHS","Con_CHS","TM_CHS"),
#     # for atrophy
#     condition_levels = c("CM","HLM","HEM","LEM","HLEM"),
#     strip_leading_token = TRUE
# ) {
#   # local helpers
#   # normalize_condition: a function to clean SAMPLE ID into a standardized condition label
#   normalize_condition_local <- function(x) {
#     x <- as.character(x)
#     base <- sub("-\\d+$", "", x, perl = TRUE)
#     if (strip_leading_token) {
#       base <- sub("^[0-9]+(?:\\.[0-9]+)?X\\s+", "", base, ignore.case = TRUE, perl = TRUE)
#       base <- sub("^[0-9]+(?:\\.[0-9]+)?X(?:_|-|\\s)+", "", base, ignore.case = TRUE, perl = TRUE)
#     }
#     base <- stringr::str_squish(base)
#     base <- gsub("[-\\s]+", "_", base, perl = TRUE)
#     up   <- toupper(base)
#     up   <- gsub("__+", "_", up, perl = TRUE)
#     dplyr::recode(up,
#                   "CM"   = "CM",
#                   "HLM"    = "HLM",
#                   "HEM" = "HEM",
#                   "LEM"  = "LEM",
#                   "HLEM" = "HLEM",
#                   # "TM_CHS"  = "TM_CHS",
#                   .default  = up)
#   }
#   # parse_ids: a function to extract condition_label and bio_replicate from SAMPLE ID
#   parse_ids_local <- function(sample_ids) {
#     rep_num <- suppressWarnings(as.integer(stringr::str_match(sample_ids, "-(\\d+)$")[,2]))
#     cond    <- normalize_condition_local(sample_ids)
#     tibble::tibble(condition_label = cond, bio_replicate = rep_num)
#   }
#   
#   # infer target genes if not provided
#   default_exclude <- c(sample_id_col, actin_col, gapdh_col, ".__sheet__")
#   tg <- if (is.null(target_gene_cols)) setdiff(names(df), default_exclude) else target_gene_cols
#   
#   # attach parsed IDs and housekeeping Cts
#   id_parsed <- parse_ids_local(df[[sample_id_col]])
#   df2 <- dplyr::bind_cols(df, id_parsed) %>%
#     dplyr::group_by(condition_label) %>%
#     dplyr::mutate(bio_replicate = ifelse(is.na(bio_replicate), dplyr::row_number(), bio_replicate)) %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(
#       ct_actin = .[[actin_col]],
#       ct_gapdh = .[[gapdh_col]]
#     )
#   
#   # pivot to long, drop missing Ct_target, tidy labels
#   long <- df2 %>%
#     dplyr::select(.__sheet__, dplyr::all_of(sample_id_col), condition_label, bio_replicate,
#                   ct_actin, ct_gapdh, dplyr::all_of(tg)) %>%
#     tidyr::pivot_longer(
#       cols = dplyr::all_of(tg),
#       names_to = "target_gene",
#       values_to = "ct_target"
#     ) %>%
#     dplyr::filter(!is.na(ct_target)) %>%
#     dplyr::mutate(condition_label = as.character(condition_label) %>% stringr::str_squish())
#   
#   # apply desired ordering (present first, others appended)
#   present <- unique(long$condition_label)
#   ordered <- c(condition_levels[condition_levels %in% present], setdiff(present, condition_levels))
#   long$condition_label <- factor(long$condition_label, levels = ordered)
#   
#   # assign technical replicate index per (gene, condition, bio) if duplicates exist
#   long <- long %>%
#     dplyr::group_by(target_gene, condition_label, bio_replicate) %>%
#     dplyr::arrange(.data[[sample_id_col]], .by_group = TRUE) %>%
#     dplyr::mutate(tech_replicate = dplyr::row_number()) %>%
#     dplyr::ungroup()
#   
#   long
# }

## build long-format Ct tables for all sheets
base_long_list <- purrr::map(
  raw_list,
  ~ build_base_long_unified(
    df = .x,
    sample_id_col = sample_id_col,
    actin_col     = actin_col,
    gapdh_col     = gapdh_col,
    target_gene_cols = target_gene_cols,
    # condition_levels = condition_levels,
    strip_leading_token = strip_leading_token))

## build_target_ct_table: a function to compute ΔCt, ΔΔCt, 2^-ΔΔCt, tests, and summaries for a chosen housekeeping gene, handling missing data/groups/calibrator per sheet and per gene
build_target_ct_table <- function(long_df, hk = c("ACTIN","GAPDH")) {
  hk <- match.arg(hk)
  
  df <- long_df %>%
    mutate(ct_housekeeping = ifelse(hk == "ACTIN", ct_actin, ct_gapdh)) %>%
    filter(!is.na(ct_housekeeping)) %>%
    mutate(delta_ct = ct_target - ct_housekeeping) %>%
    { if (!"tech_replicate" %in% names(.)) dplyr::mutate(., tech_replicate = 1L) else . }
  
  ## create your mapped gene and housekeeping data
  qpcr_replicate_data <- df %>%
    mutate(
      condition_label = as.character(condition_label) %>% str_squish()
    )
  
  condition_levels <- qpcr_replicate_data %>%
    distinct(condition_label) %>%
    pull(condition_label)
  
  qpcr_replicate_data <- qpcr_replicate_data %>%
    mutate(condition_label = factor(condition_label, levels = condition_levels))
  
  ## Define biological and technical replicates
  qpcr_replicate_data <- qpcr_replicate_data %>%
    dplyr::group_by(target_gene, condition_label) %>%
    dplyr::mutate(
      row_in_condition = dplyr::row_number(),
      n_rows_condition = dplyr::n()
      # keep existing bio_replicate and tech_replicate as-is
    ) %>%
    dplyr::ungroup()
  
  ## collapse technical duplicates to biological-level data
  qpcr_bio_data <- qpcr_replicate_data %>%
    group_by(target_gene, condition_label, bio_replicate) %>%
    summarise(
      ct_target       = mean(ct_target,       na.rm = TRUE),
      ct_housekeeping = mean(ct_housekeeping, na.rm = TRUE),
      delta_ct        = mean(delta_ct,        na.rm = TRUE),
      .groups         = "drop")
  
  ## Compute ΔΔCt and fold change
  qpcr_bio_data <- qpcr_bio_data %>%
    group_by(target_gene) %>%
    mutate(
      has_calibrator      = any(condition_label == calibrator_condition_label & !is.na(delta_ct)),
      calibrator_delta_ct = ifelse(
        has_calibrator,
        mean(delta_ct[condition_label == calibrator_condition_label], na.rm = TRUE),
        NA_real_
      )
    ) %>%
    ungroup() %>%
    mutate(
      delta_delta_ct = ifelse(is.na(calibrator_delta_ct), NA_real_, delta_ct - calibrator_delta_ct),
      fold_change    = ifelse(is.na(delta_delta_ct),      NA_real_, 2^(-delta_delta_ct))
    )
  
  ## Perform statistical testing
  gene_condition_tests <- qpcr_bio_data %>%
    group_by(target_gene, condition_label) %>%
    reframe({
      fc <- fold_change
      fc <- fc[!is.na(fc)]
      n  <- length(fc)
      
      normal <- if (n >= 3) {
        tryCatch(shapiro.test(fc)$p.value > 0.05, error = function(e) FALSE)
      } else FALSE
      
      p_val <- NA_real_
      if (n == 0) {
        p_val <- NA_real_
      } else if (normal && n >= 2) {
        p_val <- tryCatch(t.test(fc, mu = 1)$p.value, error = function(e) NA_real_)
      } else {
        all_equal_mu <- all(abs(fc - 1) < .Machine$double.eps^0.5)
        if (all_equal_mu) {
          p_val <- 1
        } else {
          p_val <- tryCatch(wilcox.test(fc, mu = 1, exact = FALSE)$p.value, error = function(e) NA_real_)
        }
      }
      
      tibble(
        normal_distribution = normal,
        p_value             = p_val,
        replicate_count     = n
      )
    })
  
  ## get mean ± SE fold change per group, merges stats results, and get sig
  gene_condition_summary <- qpcr_bio_data %>%
    group_by(target_gene, condition_label) %>%
    summarise(
      mean_fold_change = mean(fold_change, na.rm = TRUE),
      n_non_na         = sum(!is.na(fold_change)),
      se_fold_change   = ifelse(n_non_na > 0, sd(fold_change, na.rm = TRUE) / sqrt(n_non_na), NA_real_),
      .groups = "drop"
    ) %>%
    select(-n_non_na) %>%
    left_join(gene_condition_tests,
              by = c("target_gene", "condition_label")) %>%
    mutate(
      significance_symbol = case_when(
        is.na(p_value)  ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ ""),
      annotation_y_pos = mean_fold_change + se_fold_change * 1.2)
  
  list(
    qpcr_replicate_data    = qpcr_replicate_data,
    qpcr_bio_data          = qpcr_bio_data,
    gene_condition_summary = gene_condition_summary
  )
}

## run build_target_ct_table per sheet, and create at least four data frames of results
results_by_sheet <- list()

for (sh in names(base_long_list)) {
  long_df <- base_long_list[[sh]]
  res_actin <- build_target_ct_table(long_df, "ACTIN")
  res_gapdh <- build_target_ct_table(long_df, "GAPDH")
  results_by_sheet[[sh]] <- list(ACTIN = res_actin, GAPDH = res_gapdh)
  
  # expose named summary data frames per your requirement
  sh_safe <- sanitize_sheet_name(sh)
  assign(paste0("gene_condition_summary_actin_", sh_safe), res_actin$gene_condition_summary, envir = .GlobalEnv)
  assign(paste0("gene_condition_summary_gapdh_", sh_safe), res_gapdh$gene_condition_summary, envir = .GlobalEnv)
}

## Defaults for downstream plotting/export: use GAPDH from the first sheet
first_sheet <- names(results_by_sheet)[1]
qpcr_replicate_data <- results_by_sheet[[first_sheet]]$GAPDH$qpcr_replicate_data
qpcr_bio_data       <- results_by_sheet[[first_sheet]]$GAPDH$qpcr_bio_data
gene_condition_summary <- results_by_sheet[[first_sheet]]$GAPDH$gene_condition_summary

# test_1 <- results_by_sheet[[first_sheet]]$ACTIN$qpcr_replicate_data
# test_2 <- results_by_sheet[[first_sheet]]$GAPDH$qpcr_replicate_data
# test_3 <- results_by_sheet[[first_sheet]]$ACTIN$qpcr_bio_data
# test_4 <- results_by_sheet[[first_sheet]]$GAPDH$qpcr_bio_data
# test_5 <- results_by_sheet[[first_sheet]]$ACTIN$gene_condition_summary
# test_6 <- results_by_sheet[[first_sheet]]$GAPDH$gene_condition_summary

# gene_condition_summary_gapdh_2_5x
# gene_condition_summary_gapdh_5x
# gene_condition_summary_actin_2_5x
# gene_condition_summary_actin_5x

# gene_condition_summary_gapdh_2_5x
# gene_condition_summary_gapdh_5x
# gene_condition_summary_gapdh_10x
# gene_condition_summary_actin_2_5x
# gene_condition_summary_actin_5x
# gene_condition_summary_actin_10x

## Plot the results  ======================

## plot_fold_change: a function to plot fc bars with error bars and significance for one gene or all genes
plot_fold_change <- function(target_gene_name = NULL) {
  
  plotting_data <- if (is.null(target_gene_name)) {
    gene_condition_summary
  } else {
    dplyr::filter(gene_condition_summary, target_gene == target_gene_name)
  }
  
  ggplot(plotting_data, aes(x = condition_label, y = mean_fold_change)) +
    geom_col(fill = "grey30", color = "grey30", width = 0.6) +
    geom_errorbar(
      aes(ymin = mean_fold_change - se_fold_change,
          ymax = mean_fold_change + se_fold_change),
      width = 0.2,
      color = "grey30"
    ) +
    geom_text(
      aes(y = annotation_y_pos, label = significance_symbol),
      vjust = 0
    ) +
    labs(
      x     = NULL,
      y     = "Relative expression (2^-ΔΔCt)",
      title = paste0("Calibrator condition: ", calibrator_condition_label)
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# plot_fold_change("FKBP5")

ggplot(gene_condition_summary, aes(x = condition_label, y = mean_fold_change)) +
  geom_col(fill = "grey30", color = "grey30", width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_fold_change - se_fold_change,
        ymax = mean_fold_change + se_fold_change),
    width = 0.2, color = "grey30") +
  geom_text(
    aes(y = annotation_y_pos, label = significance_symbol),
    vjust = 0) +
  facet_wrap(~ target_gene, scales = "free_y") +
  labs(x = NULL, y = "Relative expression (2^-ΔΔCt)", title = paste0("Calibrator cond: ", calibrator_condition_label)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## OR do a faceted plot

## for spleen
## make sure the groups are ordered as condition_levels

gene_condition_summary_gapdh_2_5x$condition_label <- factor(gene_condition_summary_gapdh_2_5x$condition_label, levels = c(condition_levels))
# gene_condition_summary_gapdh_5x$condition_label <- factor(gene_condition_summary_gapdh_5x$condition_label, levels = c(condition_levels))
# gene_condition_summary_gapdh_10x$condition_label <- factor(gene_condition_summary_gapdh_10x$condition_label, levels = c(condition_levels))

gene_condition_summary_actin_2_5x$condition_label <- factor(gene_condition_summary_actin_2_5x$condition_label, levels = c(condition_levels))
# gene_condition_summary_actin_5x$condition_label <- factor(gene_condition_summary_actin_5x$condition_label, levels = c(condition_levels))
# gene_condition_summary_actin_10x$condition_label <- factor(gene_condition_summary_actin_10x$condition_label, levels = c(condition_levels))

# ggplot(gene_condition_summary_gapdh_2_5x, aes(x = condition_label, y = mean_fold_change)) +
#   geom_col(fill = "grey30", color = "grey30", width = 0.6) +
#   geom_errorbar(
#     aes(ymin = mean_fold_change - se_fold_change,
#         ymax = mean_fold_change + se_fold_change),
#     width = 0.2, color = "grey30") +
#   geom_text(
#     aes(y = annotation_y_pos, label = significance_symbol),
#     vjust = 0) +
#   facet_wrap(~ target_gene, scales = "free_y") +
#   labs(x = NULL, y = "Relative expression (2^-ΔΔCt)", title = paste0("Calibrator cond: ", calibrator_condition_label)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ggplot(gene_condition_summary_gapdh_5x, aes(x = condition_label, y = mean_fold_change)) +
#   geom_col(fill = "grey30", color = "grey30", width = 0.6) +
#   geom_errorbar(
#     aes(ymin = mean_fold_change - se_fold_change,
#         ymax = mean_fold_change + se_fold_change),
#     width = 0.2, color = "grey30") +
#   geom_text(
#     aes(y = annotation_y_pos, label = significance_symbol),
#     vjust = 0) +
#   facet_wrap(~ target_gene, scales = "free_y") +
#   labs(x = NULL, y = "Relative expression (2^-ΔΔCt)", title = paste0("Calibrator cond: ", calibrator_condition_label)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ggplot(gene_condition_summary_gapdh_10x, aes(x = condition_label, y = mean_fold_change)) +
#   geom_col(fill = "grey30", color = "grey30", width = 0.6) +
#   geom_errorbar(
#     aes(ymin = mean_fold_change - se_fold_change,
#         ymax = mean_fold_change + se_fold_change),
#     width = 0.2, color = "grey30") +
#   geom_text(
#     aes(y = annotation_y_pos, label = significance_symbol),
#     vjust = 0) +
#   facet_wrap(~ target_gene, scales = "free_y") +
#   labs(x = NULL, y = "Relative expression (2^-ΔΔCt)", title = paste0("Calibrator cond: ", calibrator_condition_label)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ggplot(gene_condition_summary_actin_2_5x, aes(x = condition_label, y = mean_fold_change)) +
#   geom_col(fill = "grey30", color = "grey30", width = 0.6) +
#   geom_errorbar(
#     aes(ymin = mean_fold_change - se_fold_change,
#         ymax = mean_fold_change + se_fold_change),
#     width = 0.2, color = "grey30") +
#   geom_text(
#     aes(y = annotation_y_pos, label = significance_symbol),
#     vjust = 0) +
#   facet_wrap(~ target_gene, scales = "free_y", ncol = 5) +
#   labs(x = NULL, y = "Relative expression (2^-ΔΔCt)", title = paste0("Calibrator cond: ", calibrator_condition_label)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ggplot(gene_condition_summary_actin_5x, aes(x = condition_label, y = mean_fold_change)) +
#   geom_col(fill = "grey30", color = "grey30", width = 0.6) +
#   geom_errorbar(
#     aes(ymin = mean_fold_change - se_fold_change,
#         ymax = mean_fold_change + se_fold_change),
#     width = 0.2, color = "grey30") +
#   geom_text(
#     aes(y = annotation_y_pos, label = significance_symbol),
#     vjust = 0) +
#   facet_wrap(~ target_gene, scales = "free_y") +
#   labs(x = NULL, y = "Relative expression (2^-ΔΔCt)", title = paste0("Calibrator cond: ", calibrator_condition_label)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ggplot(gene_condition_summary_actin_10x, aes(x = condition_label, y = mean_fold_change)) +
#   geom_col(fill = "grey30", color = "grey30", width = 0.6) +
#   geom_errorbar(
#     aes(ymin = mean_fold_change - se_fold_change,
#         ymax = mean_fold_change + se_fold_change),
#     width = 0.2, color = "grey30") +
#   geom_text(
#     aes(y = annotation_y_pos, label = significance_symbol),
#     vjust = 0) +
#   facet_wrap(~ target_gene, scales = "free_y") +
#   labs(x = NULL, y = "Relative expression (2^-ΔΔCt)", title = paste0("Calibrator cond: ", calibrator_condition_label)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Save per gene tables to Excel - v1  ======================

# ## v1 - per primer concentration (aka, gene_condition_summary)
# 
# ## for export I kept the original workflow and used the defaults (GAPDH, first sheet).
# ## if one wants to export for another sheet/HK, reassign qpcr_* objects from results_by_sheet.
# 
# gene_condition_summary <- results_by_sheet[[first_sheet]]$GAPDH$gene_condition_summary
# 
# ## wide table with technical duplicates per biological replicate
# ## (in wide input Ct2/CT2 may be NA)
# replicate_wide <- qpcr_replicate_data %>%
#   group_by(target_gene, condition_label, bio_replicate) %>%
#   summarise(
#     Ct1 = ct_target[tech_replicate == 1][1],
#     Ct2 = ct_target[tech_replicate == 2][1],
#     CT1 = ct_housekeeping[tech_replicate == 1][1],
#     CT2 = ct_housekeeping[tech_replicate == 2][1],
#     .groups = "drop"
#   ) %>%
#   mutate(
#     Target_gene_Ct_avg  = rowMeans(cbind(Ct1, Ct2), na.rm = TRUE),
#     Housekeeping_Ct_avg = rowMeans(cbind(CT1, CT2), na.rm = TRUE)
#   )
# 
# ## mean and SD of ΔCt and ΔΔCt per condition from biological replicates
# delta_stats_bio <- qpcr_bio_data %>%
#   group_by(target_gene, condition_label) %>%
#   summarise(
#     Mean       = mean(delta_ct,        na.rm = TRUE),
#     SD         = sd(delta_ct,          na.rm = TRUE),
#     Mean_ddCt  = mean(delta_delta_ct,  na.rm = TRUE),
#     `2^-ΔΔCt`  = ifelse(is.na(Mean_ddCt), NA_real_, 2^(-Mean_ddCt)),
#     `2^-mean`  = ifelse(is.na(Mean),     NA_real_, 2^(-Mean)),
#     .groups    = "drop")
# 
# ## build per-gene export table: biological replicates per condition, with Ct1/Ct2 and housekeeping Ct1/Ct2
# per_gene_export <- qpcr_bio_data %>%
#   left_join(replicate_wide,
#             by = c("target_gene", "condition_label", "bio_replicate")) %>%
#   left_join(delta_stats_bio,
#             by = c("target_gene", "condition_label")) %>%
#   mutate(
#     Name      = condition_label,
#     Replicate = bio_replicate,
#     `ΔCt`     = delta_ct,
#     `ΔΔCt`    = Mean_ddCt) %>%
#   select(target_gene, Name, Replicate, Ct1, Ct2, Target_gene_Ct_avg,
#          CT1, CT2, Housekeeping_Ct_avg,
#          `ΔCt`, Mean, SD, `ΔΔCt`, `2^-ΔΔCt`,`2^-mean`)
# 
# per_gene_list <- per_gene_export %>%
#   group_by(target_gene) %>%
#   group_split()
# 
# names(per_gene_list) <- per_gene_export %>%
#   distinct(target_gene) %>%
#   pull(target_gene)
# 
# per_gene_list <- lapply(per_gene_list, function(df) dplyr::select(df, -target_gene))
# 
# output_file <- sub("\\.xlsx$", "_per_gene_tables.xlsx", qpcr_file_path)
# 
# # write_xlsx(per_gene_list, path = output_file)

## Save per gene tables to Excel - v2 --------------------------------------

## v2 - for all primer concentration (aka, gene_condition_summary______)

## export all sheets × both housekeeping genes (ACTIN + GAPDH) and plots
## concentration label is auto-detected from sheet name (e.g., "2.5", "5", "10" with/without "X")

## get concentration tag from sheet name
get_conc_tag <- function(sheet_name) {
  x <- as.character(sheet_name)
  m <- stringr::str_match(x, "([0-9]+(?:\\.[0-9]+)?)\\s*[xX]?")
  ifelse(is.na(m[,2]), sanitize_sheet_name(x), gsub("\\.", "_", m[,2]))
}

## export_per_gene_tables_plots: export one excel sheet × housekeeping gene as per-gene workbook (each gene = one tab) and a matched faceted fold-change plot 
export_per_gene_tables_plots <- function(sh, hk = c("ACTIN","GAPDH")) {
  hk <- match.arg(hk)
  
  ## pull the correct objects for this sheet & HK
  qpcr_replicate_data <- results_by_sheet[[sh]][[hk]]$qpcr_replicate_data
  qpcr_bio_data       <- results_by_sheet[[sh]][[hk]]$qpcr_bio_data
  gene_condition_summary <- results_by_sheet[[sh]][[hk]]$gene_condition_summary
  
  ## wide table with technical duplicates per biological replicate
  ## (in wide input Ct2/CT2 may be NA)
  replicate_wide <- qpcr_replicate_data %>%
    group_by(target_gene, condition_label, bio_replicate) %>%
    summarise(
      Ct1 = ct_target[tech_replicate == 1][1],
      Ct2 = ct_target[tech_replicate == 2][1],
      CT1 = ct_housekeeping[tech_replicate == 1][1],
      CT2 = ct_housekeeping[tech_replicate == 2][1],
      .groups = "drop"
    ) %>%
    mutate(
      Target_gene_Ct_avg  = rowMeans(cbind(Ct1, Ct2), na.rm = TRUE),
      Housekeeping_Ct_avg = rowMeans(cbind(CT1, CT2), na.rm = TRUE)
    )
  
  ## mean and SD of ΔCt and ΔΔCt per condition from biological replicates
  delta_stats_bio <- qpcr_bio_data %>%
    group_by(target_gene, condition_label) %>%
    summarise(
      Mean_dCt            = mean(delta_ct,        na.rm = TRUE),
      SD_dct              = sd(delta_ct,          na.rm = TRUE),
      `Mean_2^-ddCt`      = mean(fold_change,     na.rm = TRUE),
      `SD_2^-ddCt`        = sd(fold_change,       na.rm = TRUE),
      .groups    = "drop")
  
  ## build per-gene export table: biological replicates per condition, with Ct1/Ct2 and housekeeping Ct1/Ct2
  per_gene_export <- qpcr_bio_data %>%
    left_join(replicate_wide,
              by = c("target_gene", "condition_label", "bio_replicate")) %>%
    left_join(delta_stats_bio,
              by = c("target_gene", "condition_label")) %>%
    mutate(
      Name      = condition_label,
      Replicate = bio_replicate,
      `ΔCt`     = delta_ct,
      `ΔΔCt`    = delta_delta_ct,
      `2^-ΔΔCt` = fold_change) %>%
    select(target_gene, Name, Replicate, Ct1, Ct2, Target_gene_Ct_avg,
           CT1, CT2, Housekeeping_Ct_avg,
           `ΔCt`, Mean_dCt, SD_dct, `ΔΔCt`,`2^-ΔΔCt`,`Mean_2^-ddCt`,`SD_2^-ddCt`)
  
  per_gene_list <- per_gene_export %>%
    group_by(target_gene) %>%
    group_split()
  
  names(per_gene_list) <- per_gene_export %>%
    distinct(target_gene) %>%
    pull(target_gene)
  
  per_gene_list <- lapply(per_gene_list, function(df) dplyr::select(df, -target_gene))
  
  ## output name: auto conc tag + HK + sheet-safe
  conc_tag <- get_conc_tag(sh)
  sh_safe  <- sanitize_sheet_name(sh)
  
  base_out <- sub("\\.xlsx$", "", qpcr_file_path)
  output_file <- paste0(base_out, "_", hk, "_", conc_tag, "x_per_gene_tables.xlsx")
  
  write_xlsx(per_gene_list, path = output_file)
  
  ## build and save matched faceted plot
  p <- ggplot(gene_condition_summary, aes(x = condition_label, y = mean_fold_change)) +
    geom_col(fill = "grey30", color = "grey30", width = 0.6) +
    geom_errorbar(
      aes(ymin = mean_fold_change - se_fold_change,
          ymax = mean_fold_change + se_fold_change),
      width = 0.2, color = "grey30") +
    geom_text(
      aes(y = annotation_y_pos, label = significance_symbol),
      vjust = 0) +
    facet_wrap(~ target_gene, scales = "free_y") +
    labs(x = NULL, y = "Relative expression (2^-ddCt)",
         title = paste0("Calibrator cond: ", calibrator_condition_label, " | ", hk, " | ", conc_tag, "x")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  plot_file <- paste0(base_out, "_", hk, "_", conc_tag, "x_per_gene_plots.pdf")
  ggsave(plot_file, plot = p, width = 10, height = 8)
  
  ## return useful objects if needed
  list(
    sheet = sh,
    hk    = hk,
    conc_tag = conc_tag,
    gene_condition_summary = gene_condition_summary,
    output_file = output_file,
    plot_file   = plot_file
  )
}

## run export_per_gene_tables_plots

export_results <- list()

for (sh in names(results_by_sheet)) {
  export_results[[paste0(sanitize_sheet_name(sh), "_ACTIN")]] <- export_per_gene_tables_plots(sh, "ACTIN")
  export_results[[paste0(sanitize_sheet_name(sh), "_GAPDH")]] <- export_per_gene_tables_plots(sh, "GAPDH")
}

## xlsx files for tables and pdf files for plots for these variables:
# gene_condition_summary_gapdh_2_5x
# gene_condition_summary_gapdh_5x
# gene_condition_summary_gapdh_10x
# gene_condition_summary_actin_2_5x
# gene_condition_summary_actin_5x
# gene_condition_summary_actin_10x

## Compare with RNA-seq results - TM_chicken_muscle [in progress] ======================

## read DEG RNA-seq results

DEG_muscle_D22_qPCR <- read_excel("C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_muscle\\TMomics_DEG_muscle_D22_qPCR_v2.xlsx", 
                                  sheet = "Sheet1")

unique(DEG_muscle_D22_qPCR$contrast)
unique(DEG_muscle_D22_qPCR$SYMBOL)
nrow(DEG_muscle_D22_qPCR)

qPCR_comp <- c("TM_N vs Con_N", "Con_AHS vs Con_N", "TM_AHS vs Con_N", "Con_CHS vs Con_N", "TM_CHS vs Con_N")

DEG_muscle_D22_qPCR_mod <- DEG_muscle_D22_qPCR %>%
  # filter(contrast %in% qPCR_comp) %>%
  select(c("SYMBOL","log2FoldChange","direction","contrast"))

unique(DEG_muscle_D22_qPCR_mod$SYMBOL)
nrow(DEG_muscle_D22_qPCR_mod)

qPCR_genes <- c("HSPH1","FKBP5","MSMB","FAM19A4","LECT2","CRB2","CBLN2","SFRP4","HSP90AA1","COL1A1","TMEM163","PLB1")
length(unique(DEG_muscle_D22_qPCR_mod$SYMBOL))
length(qPCR_genes)

## read qPCR processed results

base_dir <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_muscle\\"
read_per_gene_file <- function(path) {
  shs <- readxl::excel_sheets(path)
  dplyr::bind_rows(lapply(shs, function(s) {
    df <- readxl::read_excel(path, sheet = s)
    df$SYMBOL <- s
    df
  }))
}
f_actin_2_5x <- file.path(base_dir, "MUSCLE_CT1_CT2_cDNA_1_400_5_2_5_mod_SS_ACTIN_2_5x_per_gene_tables.xlsx")
f_actin_5x   <- file.path(base_dir, "MUSCLE_CT1_CT2_cDNA_1_400_5_2_5_mod_SS_ACTIN_5x_per_gene_tables.xlsx")
f_gapdh_2_5x <- file.path(base_dir, "MUSCLE_CT1_CT2_cDNA_1_400_5_2_5_mod_SS_GAPDH_2_5x_per_gene_tables.xlsx")
f_gapdh_5x   <- file.path(base_dir, "MUSCLE_CT1_CT2_cDNA_1_400_5_2_5_mod_SS_GAPDH_5x_per_gene_tables.xlsx")
per_gene_tables_actin_2_5x <- read_per_gene_file(f_actin_2_5x)
per_gene_tables_actin_5x   <- read_per_gene_file(f_actin_5x)
per_gene_tables_gapdh_2_5x <- read_per_gene_file(f_gapdh_2_5x)
per_gene_tables_gapdh_5x   <- read_per_gene_file(f_gapdh_5x)

sum_per_gene_tables_actin_2_5x <- summarise(group_by(per_gene_tables_actin_2_5x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")
sum_per_gene_tables_actin_5x   <- summarise(group_by(per_gene_tables_actin_5x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")
sum_per_gene_tables_gapdh_2_5x <- summarise(group_by(per_gene_tables_gapdh_2_5x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")
sum_per_gene_tables_gapdh_5x   <- summarise(group_by(per_gene_tables_gapdh_5x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")

cond_levels <- c("Con_N", "TM_N", "Con_AHS", "TM_AHS", "Con_CHS", "TM_CHS")

sum_per_gene_tables_actin_2_5x$Name <- factor(sum_per_gene_tables_actin_2_5x$Name,levels=cond_levels)
sum_per_gene_tables_actin_5x$Name <- factor(sum_per_gene_tables_actin_5x$Name,levels=cond_levels)
sum_per_gene_tables_gapdh_2_5x$Name <- factor(sum_per_gene_tables_gapdh_2_5x$Name,levels=cond_levels)
sum_per_gene_tables_gapdh_5x$Name <- factor(sum_per_gene_tables_gapdh_5x$Name,levels=cond_levels)
sum_per_gene_tables_actin_2_5x <- sum_per_gene_tables_actin_2_5x[order(sum_per_gene_tables_actin_2_5x$Name), ]
sum_per_gene_tables_actin_5x <- sum_per_gene_tables_actin_5x[order(sum_per_gene_tables_actin_5x$Name), ]
sum_per_gene_tables_gapdh_2_5x <- sum_per_gene_tables_gapdh_2_5x[order(sum_per_gene_tables_gapdh_2_5x$Name), ]
sum_per_gene_tables_gapdh_5x <- sum_per_gene_tables_gapdh_5x[order(sum_per_gene_tables_gapdh_5x$Name), ]

## go through each gene_condition_summary_ or sum_per_gene_tables_ to compare expression pattern with DEG_muscle_D22_qPCR_mod

## (1) gene_condition_summary_actin_2_5x    OR    sum_per_gene_tables_actin_2_5x
## (2) gene_condition_summary_actin_5x      OR    sum_per_gene_tables_actin_5x
## (3) gene_condition_summary_gapdh_2_5x    OR    sum_per_gene_tables_gapdh_2_5x
## (4) gene_condition_summary_gapdh_5x      OR    sum_per_gene_tables_gapdh_5x

target_qPCR_gene <- c("HSPH1")

DEG_muscle_D22_qPCR_mod[DEG_muscle_D22_qPCR_mod$SYMBOL == target_qPCR_gene,]
sum_per_gene_tables_actin_2_5x[sum_per_gene_tables_actin_2_5x$SYMBOL == target_qPCR_gene,]
sum_per_gene_tables_actin_5x[sum_per_gene_tables_actin_5x$SYMBOL == target_qPCR_gene,]

DEG_muscle_D22_qPCR_mod[DEG_muscle_D22_qPCR_mod$SYMBOL == target_qPCR_gene,]
sum_per_gene_tables_gapdh_2_5x[sum_per_gene_tables_gapdh_2_5x$SYMBOL == target_qPCR_gene,]
sum_per_gene_tables_gapdh_5x[sum_per_gene_tables_gapdh_5x$SYMBOL == target_qPCR_gene,]

## filter the matched genes/condition in gene_condition_summary_ to get the final plot

passed_qPCR_genes <- c("HSPH1","FKBP5","MSMB","LECT2","SFRP4","HSP90AA1","COL1A1")
length(passed_qPCR_genes)

passed_qPCR_genes_actin_2_5x <- c("HSPH1","FKBP5","LECT2","SFRP4","HSP90AA1")
passed_qPCR_genes_actin_5x <- c("MSMB","COL1A1")

passed_gene_condition_summary_actin_2_5x <- gene_condition_summary_actin_2_5x %>%
  filter(target_gene %in% passed_qPCR_genes_actin_2_5x)
passed_gene_condition_summary_actin_5x <- gene_condition_summary_actin_5x %>%
  filter(target_gene %in% passed_qPCR_genes_actin_5x)

passed_gene_condition_summary_actin <- rbind(passed_gene_condition_summary_actin_2_5x,
                                             passed_gene_condition_summary_actin_5x)

## plot the final fig

## plot parameters

grey_map <- c(
  Con_N   = "black",
  TM_N    = "grey25",
  Con_AHS = "grey40",
  TM_AHS  = "grey55", 
  Con_CHS = "grey45",
  TM_CHS  = "grey65")

passed_qPCR_genes_levels <- c("COL1A1","HSP90AA1",
                              "HSPH1","FKBP5",
                              "MSMB","LECT2","SFRP4")
passed_gene_condition_summary_actin$target_gene <- factor(passed_gene_condition_summary_actin$target_gene,
                                                          levels=passed_qPCR_genes_levels)
passed_gene_condition_summary_actin <- passed_gene_condition_summary_actin[order(passed_gene_condition_summary_actin$target_gene), ]

passed_qPCR_genes_levels <- list(
  set1 = c("HSPH1","FKBP5"),
  set2 = c("HSP90AA1","COL1A1"),
  set3 = c("MSMB","LECT2","SFRP4"))

## final plot

ggplot(subset(passed_gene_condition_summary_actin, 
              # target_gene %in% passed_qPCR_genes_levels$set3
              ),
       aes(x = condition_label, y = mean_fold_change,
           fill = condition_label, color = condition_label)) +
  geom_col(width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_fold_change - se_fold_change,
        ymax = mean_fold_change + se_fold_change),
    width = 0.2, color = "grey30") +
  geom_text(
    aes(y = annotation_y_pos, label = significance_symbol),
    vjust = 0) +
  scale_fill_manual(values = grey_map) +
  scale_color_manual(values = grey_map) +
  facet_wrap(~ target_gene, scales = "free", ncol = 4) +
  labs(x = NULL, y = "Relative expression (2^-ΔΔCt)") +
  theme_classic() +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none", 
    strip.background = element_blank(),
    strip.text = element_text(face = "bold",size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1))

## compare actin and gapdh HK genes results

passed_qPCR_genes_actin_2_5x <- c("HSPH1","FKBP5","LECT2","SFRP4","HSP90AA1")
passed_qPCR_genes_actin_5x <- c("MSMB","COL1A1")

sum_per_gene_tables_actin_2_5x_mod <- sum_per_gene_tables_actin_2_5x %>% filter(SYMBOL %in% passed_qPCR_genes_actin_2_5x)
sum_per_gene_tables_gapdh_2_5x_mod <- sum_per_gene_tables_gapdh_2_5x %>% filter(SYMBOL %in% passed_qPCR_genes_actin_2_5x)
sum_per_gene_tables_actin_5x_mod <- sum_per_gene_tables_actin_5x %>% filter(SYMBOL %in% passed_qPCR_genes_actin_5x)
sum_per_gene_tables_gapdh_5x_mod <- sum_per_gene_tables_gapdh_5x %>% filter(SYMBOL %in% passed_qPCR_genes_actin_5x)

nrow(sum_per_gene_tables_actin_2_5x)
nrow(sum_per_gene_tables_gapdh_2_5x)
nrow(sum_per_gene_tables_actin_5x)
nrow(sum_per_gene_tables_gapdh_5x)

nrow(sum_per_gene_tables_actin_2_5x_mod)
nrow(sum_per_gene_tables_gapdh_2_5x_mod)
nrow(sum_per_gene_tables_actin_5x_mod)
nrow(sum_per_gene_tables_gapdh_5x_mod)

sum_per_gene_tables_actin_2_5x_mod$housekeeping <- "ACTB"
sum_per_gene_tables_gapdh_2_5x_mod$housekeeping <- "GAPDH"
sum_per_gene_tables_actin_5x_mod$housekeeping <- "ACTB"
sum_per_gene_tables_gapdh_5x_mod$housekeeping <- "GAPDH"

actin_gapdh_comp_2_5x <- rbind(sum_per_gene_tables_actin_2_5x_mod, sum_per_gene_tables_gapdh_2_5x_mod)
actin_gapdh_comp_5x <- rbind(sum_per_gene_tables_actin_5x_mod, sum_per_gene_tables_gapdh_5x_mod)

hk_comp_2_5x <- actin_gapdh_comp_2_5x %>%
  select(SYMBOL, Name, fc, housekeeping) %>%
  tidyr::pivot_wider(names_from = housekeeping, values_from = fc)
hk_comp_5x <- actin_gapdh_comp_5x %>%
  select(SYMBOL, Name, fc, housekeeping) %>%
  tidyr::pivot_wider(names_from = housekeeping, values_from = fc)

hk_comp_2_5x$diff <- round((hk_comp_2_5x$ACTB - hk_comp_2_5x$GAPDH),7)
hk_comp_5x$diff <- round((hk_comp_5x$ACTB - hk_comp_5x$GAPDH),7)

summary(hk_comp_2_5x$diff)
summary(hk_comp_5x$diff)

cor(hk_comp_2_5x$ACTB, hk_comp_2_5x$GAPDH, use = "complete.obs")
cor(hk_comp_5x$ACTB, hk_comp_5x$GAPDH, use = "complete.obs")

## plot RNAseq and RT-qPCR res together

passed_qPCR_genes

DEG_muscle_D22_qPCR_mod_passed <- DEG_muscle_D22_qPCR_mod[DEG_muscle_D22_qPCR_mod$SYMBOL %in% passed_qPCR_genes,]

DEG_muscle_D22_qPCR_mod_passed_mod <- DEG_muscle_D22_qPCR_mod_passed %>%
  mutate(method = "RNA_seq") %>%
  rename(target_gene = SYMBOL,
         log2_fold_change = log2FoldChange) %>%
  mutate(condition_label = gsub(" ", "", sub(" .*", "", contrast))) %>%
  as.data.frame()

passed_gene_condition_summary_actin_mod <- passed_gene_condition_summary_actin %>%
  mutate(method = "RT_qPCR") %>%
  select(target_gene, mean_fold_change, method, condition_label, significance_symbol) %>%
  # rename(fold_change = mean_fold_change) %>%
  mutate(log2_fold_change = log2(mean_fold_change)) %>%
  as.data.frame() 
  
RNAseq_RTqPCR_df <- bind_rows(DEG_muscle_D22_qPCR_mod_passed_mod, passed_gene_condition_summary_actin_mod) %>%
  filter(condition_label != "Con_N")

ggplot(RNAseq_RTqPCR_df, aes(x = condition_label, y = log2_fold_change, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = significance_symbol), position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("RNA_seq" = "black", "RT_qPCR" = "gray")) +
  facet_wrap(~ target_gene, scales = "free", ncol = 4) +
  labs(x = NULL, y = "log2 FC") +
  theme_classic() +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold",size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1))

## Compare with RNA-seq results - TM_chicken_liver  [in progress] ======================

## read DEG RNA-seq results

DEG_liver_D22_qPCR <- read_excel("C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_liver\\TMomics_DEG_liver_D22_qPCR_v2_final.xlsx", 
                                  sheet = "Sheet1")

unique(DEG_liver_D22_qPCR$contrast)
unique(DEG_liver_D22_qPCR$SYMBOL)
nrow(DEG_liver_D22_qPCR)

qPCR_comp <- c("TM_N vs Con_N", "Con_AHS vs Con_N", "TM_AHS vs Con_N")

DEG_liver_D22_qPCR_mod <- DEG_liver_D22_qPCR %>%
  # filter(contrast %in% qPCR_comp) %>%
  select(c("SYMBOL","log2FoldChange","direction","contrast"))

unique(DEG_liver_D22_qPCR_mod$SYMBOL)
nrow(DEG_liver_D22_qPCR_mod)

qPCR_genes <- c("SPHKAP", "HSPH1", "DNAJA4", "COL15A1", "LOC107049986", "GADD45G", "VGLL3")
length(unique(DEG_liver_D22_qPCR_mod$SYMBOL))
length(qPCR_genes)

## read qPCR processed results

base_dir <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_liver\\"
read_per_gene_file <- function(path) {
  shs <- readxl::excel_sheets(path)
  dplyr::bind_rows(lapply(shs, function(s) {
    df <- readxl::read_excel(path, sheet = s)
    df$SYMBOL <- s
    df
  }))
}
f_actin_2_5x <- file.path(base_dir, "LIVER_CT1_CT2_cDNA_1_XXX_5_2_5_SS_ACTIN_2_5x_per_gene_tables.xlsx")
f_actin_5x   <- file.path(base_dir, "LIVER_CT1_CT2_cDNA_1_XXX_5_2_5_SS_ACTIN_5x_per_gene_tables.xlsx")
f_gapdh_2_5x <- file.path(base_dir, "LIVER_CT1_CT2_cDNA_1_XXX_5_2_5_SS_GAPDH_2_5x_per_gene_tables.xlsx")
f_gapdh_5x   <- file.path(base_dir, "LIVER_CT1_CT2_cDNA_1_XXX_5_2_5_SS_GAPDH_5x_per_gene_tables.xlsx")
per_gene_tables_actin_2_5x <- read_per_gene_file(f_actin_2_5x)
per_gene_tables_actin_5x   <- read_per_gene_file(f_actin_5x)
per_gene_tables_gapdh_2_5x <- read_per_gene_file(f_gapdh_2_5x)
per_gene_tables_gapdh_5x   <- read_per_gene_file(f_gapdh_5x)

sum_per_gene_tables_actin_2_5x <- summarise(group_by(per_gene_tables_actin_2_5x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")
sum_per_gene_tables_actin_5x   <- summarise(group_by(per_gene_tables_actin_5x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")
sum_per_gene_tables_gapdh_2_5x <- summarise(group_by(per_gene_tables_gapdh_2_5x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")
sum_per_gene_tables_gapdh_5x   <- summarise(group_by(per_gene_tables_gapdh_5x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")

cond_levels <- c("Con_N", "TM_N", "Con_AHS", "TM_AHS")

sum_per_gene_tables_actin_2_5x$Name <- factor(sum_per_gene_tables_actin_2_5x$Name,levels=cond_levels)
sum_per_gene_tables_actin_5x$Name <- factor(sum_per_gene_tables_actin_5x$Name,levels=cond_levels)
sum_per_gene_tables_gapdh_2_5x$Name <- factor(sum_per_gene_tables_gapdh_2_5x$Name,levels=cond_levels)
sum_per_gene_tables_gapdh_5x$Name <- factor(sum_per_gene_tables_gapdh_5x$Name,levels=cond_levels)
sum_per_gene_tables_actin_2_5x <- sum_per_gene_tables_actin_2_5x[order(sum_per_gene_tables_actin_2_5x$Name), ]
sum_per_gene_tables_actin_5x <- sum_per_gene_tables_actin_5x[order(sum_per_gene_tables_actin_5x$Name), ]
sum_per_gene_tables_gapdh_2_5x <- sum_per_gene_tables_gapdh_2_5x[order(sum_per_gene_tables_gapdh_2_5x$Name), ]
sum_per_gene_tables_gapdh_5x <- sum_per_gene_tables_gapdh_5x[order(sum_per_gene_tables_gapdh_5x$Name), ]

## go through each gene_condition_summary_ or sum_per_gene_tables_ to compare expression pattern with DEG_liver_D22_qPCR_mod

## (1) gene_condition_summary_actin_2_5x    OR    sum_per_gene_tables_actin_2_5x
## (2) gene_condition_summary_actin_5x      OR    sum_per_gene_tables_actin_5x
## (3) gene_condition_summary_gapdh_2_5x    OR    sum_per_gene_tables_gapdh_2_5x
## (4) gene_condition_summary_gapdh_5x      OR    sum_per_gene_tables_gapdh_5x

qPCR_genes <- c("SPHKAP", "HSPH1", "DNAJA4", "COL15A1", "LOC107049986", "GADD45G", "VGLL3")

target_qPCR_gene <- c("GADD45G")

DEG_liver_D22_qPCR_mod[DEG_liver_D22_qPCR_mod$SYMBOL == target_qPCR_gene,]
sum_per_gene_tables_actin_2_5x[sum_per_gene_tables_actin_2_5x$SYMBOL == target_qPCR_gene,]
sum_per_gene_tables_actin_5x[sum_per_gene_tables_actin_5x$SYMBOL == target_qPCR_gene,]

DEG_liver_D22_qPCR_mod[DEG_liver_D22_qPCR_mod$SYMBOL == target_qPCR_gene,]
sum_per_gene_tables_gapdh_2_5x[sum_per_gene_tables_gapdh_2_5x$SYMBOL == target_qPCR_gene,]
sum_per_gene_tables_gapdh_5x[sum_per_gene_tables_gapdh_5x$SYMBOL == target_qPCR_gene,]

## filter the matched genes/condition in gene_condition_summary_ to get the final plot

passed_qPCR_genes <- c("HSPH1", "DNAJA4", "LOC107049986", "GADD45G", "COL15A1")

length(passed_qPCR_genes)

passed_qPCR_genes_actin_2_5x <- c("HSPH1","DNAJA4","LOC107049986","GADD45G")
passed_qPCR_genes_actin_5x <- c("COL15A1")

passed_gene_condition_summary_actin_2_5x <- gene_condition_summary_actin_2_5x %>%
  filter(target_gene %in% passed_qPCR_genes_actin_2_5x)
passed_gene_condition_summary_actin_5x <- gene_condition_summary_actin_5x %>%
  filter(target_gene %in% passed_qPCR_genes_actin_5x)

passed_gene_condition_summary_actin <- rbind(passed_gene_condition_summary_actin_2_5x,
                                             passed_gene_condition_summary_actin_5x)

## the case of GADD45G
## GADD45G has an issue as it's not been run for Con_N, which means has no calibrator producing NA
## but actually it's consistent with TM_AHS and Con_AHS comp, so we will report via excel calcu.

passed_gene_condition_summary_actin <- passed_gene_condition_summary_actin %>%
  mutate(
    mean_fold_change = case_when(
      target_gene == "GADD45G" & condition_label == "Con_AHS" ~ 1.039928289,
      target_gene == "GADD45G" & condition_label == "TM_AHS"  ~ 1.296372591,
      TRUE ~ mean_fold_change),
    se_fold_change = case_when(
      target_gene == "GADD45G" & condition_label == "Con_AHS" ~ 0.188959345,
      target_gene == "GADD45G" & condition_label == "TM_AHS"  ~ 0.140672956,
      TRUE ~ se_fold_change))

## plot the final fig

## plot parameters

grey_map <- c(
  Con_N   = "black",
  TM_N    = "grey30",
  Con_AHS = "grey55",
  TM_AHS  = "grey70")

passed_gene_condition_summary_actin <- passed_gene_condition_summary_actin %>%
  mutate(target_gene = paste0("italic('", target_gene, "')"))

passed_qPCR_genes_levels <- c("italic('HSPH1')",
                              "italic('DNAJA4')",
                              "italic('LOC107049986')",
                              "italic('COL15A1')",
                              "italic('GADD45G')")

passed_gene_condition_summary_actin$target_gene <- factor(passed_gene_condition_summary_actin$target_gene,
                                                          levels=passed_qPCR_genes_levels)
passed_gene_condition_summary_actin <- passed_gene_condition_summary_actin[order(passed_gene_condition_summary_actin$target_gene), ]

passed_qPCR_genes_levels <- list(
  set1 = c("italic('HSPH1')", "italic('DNAJA4')", "italic('LOC107049986')", "italic('COL15A1')"),
  set2 = c("italic('GADD45G')"))

## final plot

ggplot(subset(passed_gene_condition_summary_actin, 
              target_gene %in% passed_qPCR_genes_levels$set1
              ),
aes(x = condition_label, y = mean_fold_change,
    fill = condition_label, color = condition_label)) +
  geom_col(width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_fold_change - se_fold_change,
        ymax = mean_fold_change + se_fold_change),
    width = 0.2, color = "black") +
  geom_text(
    aes(y = annotation_y_pos, label = significance_symbol),
    vjust = 0, color = "black") +
  scale_fill_manual(values = grey_map) +
  scale_color_manual(values = grey_map) +
  facet_wrap(~ target_gene, scales = "free", ncol = 5, labeller = label_parsed) +
  labs(x = NULL, y = "Relative expression (2^-ΔΔCt)") +
  theme_classic() +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none", 
    strip.background = element_blank(),
    strip.text = element_text(face = "bold",size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1))

## compare actin and gapdh HK genes results

passed_qPCR_genes_actin_2_5x <- c("HSPH1","DNAJA4","LOC107049986","GADD45G")
passed_qPCR_genes_actin_5x <- c("COL15A1")

sum_per_gene_tables_actin_2_5x_mod <- sum_per_gene_tables_actin_2_5x %>% filter(SYMBOL %in% passed_qPCR_genes_actin_2_5x)
sum_per_gene_tables_gapdh_2_5x_mod <- sum_per_gene_tables_gapdh_2_5x %>% filter(SYMBOL %in% passed_qPCR_genes_actin_2_5x)
sum_per_gene_tables_actin_5x_mod <- sum_per_gene_tables_actin_5x %>% filter(SYMBOL %in% passed_qPCR_genes_actin_5x)
sum_per_gene_tables_gapdh_5x_mod <- sum_per_gene_tables_gapdh_5x %>% filter(SYMBOL %in% passed_qPCR_genes_actin_5x)

nrow(sum_per_gene_tables_actin_2_5x)
nrow(sum_per_gene_tables_gapdh_2_5x)
nrow(sum_per_gene_tables_actin_5x)
nrow(sum_per_gene_tables_gapdh_5x)

nrow(sum_per_gene_tables_actin_2_5x_mod)
nrow(sum_per_gene_tables_gapdh_2_5x_mod)
nrow(sum_per_gene_tables_actin_5x_mod)
nrow(sum_per_gene_tables_gapdh_5x_mod)

sum_per_gene_tables_actin_2_5x_mod$housekeeping <- "ACTB"
sum_per_gene_tables_gapdh_2_5x_mod$housekeeping <- "GAPDH"
sum_per_gene_tables_actin_5x_mod$housekeeping <- "ACTB"
sum_per_gene_tables_gapdh_5x_mod$housekeeping <- "GAPDH"

actin_gapdh_comp_2_5x <- rbind(sum_per_gene_tables_actin_2_5x_mod, sum_per_gene_tables_gapdh_2_5x_mod)
actin_gapdh_comp_5x <- rbind(sum_per_gene_tables_actin_5x_mod, sum_per_gene_tables_gapdh_5x_mod)

hk_comp_2_5x <- actin_gapdh_comp_2_5x %>%
  select(SYMBOL, Name, fc, housekeeping) %>%
  tidyr::pivot_wider(names_from = housekeeping, values_from = fc)
hk_comp_5x <- actin_gapdh_comp_5x %>%
  select(SYMBOL, Name, fc, housekeeping) %>%
  tidyr::pivot_wider(names_from = housekeeping, values_from = fc)

hk_comp_2_5x$diff <- round((hk_comp_2_5x$ACTB - hk_comp_2_5x$GAPDH),7)
hk_comp_5x$diff <- round((hk_comp_5x$ACTB - hk_comp_5x$GAPDH),7)

summary(hk_comp_2_5x$diff)
summary(hk_comp_5x$diff)

cor(hk_comp_2_5x$ACTB, hk_comp_2_5x$GAPDH, use = "complete.obs")
cor(hk_comp_5x$ACTB, hk_comp_5x$GAPDH, use = "complete.obs")

plot(log2(hk_comp_2_5x$ACTB),log2(hk_comp_2_5x$GAPDH), pch = 16)
abline(lm(log2(hk_comp_2_5x$GAPDH) ~ log2(hk_comp_2_5x$ACTB)), col = "firebrick", lwd = 2)

## quick and dirty model of X vs 2^-X

x <- seq(-50, 50, by = 0.5)
y <- 2^(-x)

model_expo <- data.frame(X = x, Y = y)

ggplot(model_expo, aes(x = X, y = Y)) +
  geom_line(size = 1) +
  labs(x = "X", y = expression(2^{-X})) +
  theme_classic()

## Compare with RNA-seq results - TM_chicken_spleen [in progress] ======================

## read DEG RNA-seq results

# DEG_spleen_D22_qPCR <- read_excel("C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_spleen\\TMomics_DEG_spleen_D22_qPCR_v2_final.xlsx", 
#                                  sheet = "Sheet1")

DEG_spleen_D22_qPCR <- read_excel("C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_spleen\\TMomics_DEG_spleen_D22_allanno_allcomps_v6_final.xlsx", 
                                  sheet = "Sheet1")

unique(DEG_spleen_D22_qPCR$contrast)
unique(DEG_spleen_D22_qPCR$SYMBOL)
nrow(DEG_spleen_D22_qPCR)

qPCR_comp <- c("TM_N vs Con_N", "Con_AHS vs Con_N", "TM_AHS vs Con_N")

DEG_spleen_D22_qPCR_mod <- DEG_spleen_D22_qPCR %>%
  # filter(contrast %in% qPCR_comp) %>%
  select(c("SYMBOL","log2FoldChange","direction","contrast"))

unique(DEG_spleen_D22_qPCR_mod$SYMBOL)
nrow(DEG_spleen_D22_qPCR_mod)
unique(DEG_spleen_D22_qPCR_mod$contrast)

qPCR_genes <- c("HSPH1","DNAJA4","NEDD4L","ZBTB32","PCK1","CKMT2","ACSBG1","SPIN1L","FKBP5","UBAP2L2","ASB1","RAD51C")
length(unique(DEG_spleen_D22_qPCR_mod$SYMBOL))
length(qPCR_genes)

## read qPCR processed results

base_dir <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_spleen\\"
read_per_gene_file <- function(path) {
  shs <- readxl::excel_sheets(path)
  dplyr::bind_rows(lapply(shs, function(s) {
    df <- readxl::read_excel(path, sheet = s)
    df$SYMBOL <- s
    df
  }))
}

## 1/200
f_actin_2_5x <- file.path(base_dir, "CT1_CT2_cDNA_1_200_mod_SS_ACTIN_2_5x_per_gene_tables.xlsx")
f_actin_5x   <- file.path(base_dir, "CT1_CT2_cDNA_1_200_mod_SS_ACTIN_5x_per_gene_tables.xlsx")
f_actin_10x  <- file.path(base_dir, "CT1_CT2_cDNA_1_200_mod_SS_ACTIN_10x_per_gene_tables.xlsx")
f_gapdh_2_5x <- file.path(base_dir, "CT1_CT2_cDNA_1_200_mod_SS_GAPDH_2_5x_per_gene_tables.xlsx")
f_gapdh_5x   <- file.path(base_dir, "CT1_CT2_cDNA_1_200_mod_SS_GAPDH_5x_per_gene_tables.xlsx")
f_gapdh_10x  <- file.path(base_dir, "CT1_CT2_cDNA_1_200_mod_SS_GAPDH_10x_per_gene_tables.xlsx")
## 1/400
f_actin_2_5x <- file.path(base_dir, "CT1_CT2_cDNA_1_400_mod_SS_ACTIN_2_5x_per_gene_tables.xlsx")
f_actin_5x   <- file.path(base_dir, "CT1_CT2_cDNA_1_400_mod_SS_ACTIN_5x_per_gene_tables.xlsx")
f_gapdh_2_5x <- file.path(base_dir, "CT1_CT2_cDNA_1_400_mod_SS_GAPDH_2_5x_per_gene_tables.xlsx")
f_gapdh_5x   <- file.path(base_dir, "CT1_CT2_cDNA_1_400_mod_SS_GAPDH_5x_per_gene_tables.xlsx")

per_gene_tables_actin_2_5x <- read_per_gene_file(f_actin_2_5x)
per_gene_tables_actin_5x   <- read_per_gene_file(f_actin_5x)
per_gene_tables_actin_10x  <- read_per_gene_file(f_actin_10x)

per_gene_tables_gapdh_2_5x <- read_per_gene_file(f_gapdh_2_5x)
per_gene_tables_gapdh_5x   <- read_per_gene_file(f_gapdh_5x)
per_gene_tables_gapdh_10x  <- read_per_gene_file(f_gapdh_10x)

sum_per_gene_tables_actin_2_5x <- summarise(group_by(per_gene_tables_actin_2_5x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")
sum_per_gene_tables_actin_5x   <- summarise(group_by(per_gene_tables_actin_5x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")
sum_per_gene_tables_actin_10x  <- summarise(group_by(per_gene_tables_actin_10x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")

sum_per_gene_tables_gapdh_2_5x <- summarise(group_by(per_gene_tables_gapdh_2_5x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")
sum_per_gene_tables_gapdh_5x   <- summarise(group_by(per_gene_tables_gapdh_5x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")
sum_per_gene_tables_gapdh_10x  <- summarise(group_by(per_gene_tables_gapdh_10x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")

cond_levels <- c("Con_N", "TM_N", "Con_AHS", "TM_AHS", "Con_CHS", "TM_CHS")

sum_per_gene_tables_actin_2_5x$Name <- factor(sum_per_gene_tables_actin_2_5x$Name,levels=cond_levels)
sum_per_gene_tables_actin_5x$Name <- factor(sum_per_gene_tables_actin_5x$Name,levels=cond_levels)
sum_per_gene_tables_actin_10x$Name <- factor(sum_per_gene_tables_actin_10x$Name,levels=cond_levels)

sum_per_gene_tables_gapdh_2_5x$Name <- factor(sum_per_gene_tables_gapdh_2_5x$Name,levels=cond_levels)
sum_per_gene_tables_gapdh_5x$Name <- factor(sum_per_gene_tables_gapdh_5x$Name,levels=cond_levels)
sum_per_gene_tables_gapdh_10x$Name <- factor(sum_per_gene_tables_gapdh_10x$Name,levels=cond_levels)

sum_per_gene_tables_actin_2_5x <- sum_per_gene_tables_actin_2_5x[order(sum_per_gene_tables_actin_2_5x$Name), ]
sum_per_gene_tables_actin_5x <- sum_per_gene_tables_actin_5x[order(sum_per_gene_tables_actin_5x$Name), ]
sum_per_gene_tables_actin_10x <- sum_per_gene_tables_actin_10x[order(sum_per_gene_tables_actin_10x$Name), ]

sum_per_gene_tables_gapdh_2_5x <- sum_per_gene_tables_gapdh_2_5x[order(sum_per_gene_tables_gapdh_2_5x$Name), ]
sum_per_gene_tables_gapdh_5x <- sum_per_gene_tables_gapdh_5x[order(sum_per_gene_tables_gapdh_5x$Name), ]
sum_per_gene_tables_gapdh_10x <- sum_per_gene_tables_gapdh_10x[order(sum_per_gene_tables_gapdh_10x$Name), ]

## go through each gene_condition_summary_ or sum_per_gene_tables_ to compare expression pattern with DEG_spleen_D22_qPCR_mod

## (1) gene_condition_summary_actin_2_5x    OR    sum_per_gene_tables_actin_2_5x
## (2) gene_condition_summary_actin_5x      OR    sum_per_gene_tables_actin_5x
## (3) gene_condition_summary_actin_10x     OR    sum_per_gene_tables_actin_10x

## (4) gene_condition_summary_gapdh_2_5x    OR    sum_per_gene_tables_gapdh_2_5x
## (5) gene_condition_summary_gapdh_5x      OR    sum_per_gene_tables_gapdh_5x
## (6) gene_condition_summary_gapdh_10x     OR    sum_per_gene_tables_gapdh_10x

qPCR_genes <- c("HSPH1","DNAJA4",
                "NEDD4L","ZBTB32",
                "PCK1","CKMT2","ACSBG1","SPIN1L",
                "FKBP5","UBAP2L2","ASB1","RAD51C")

target_qPCR_gene <- c("PCK1") 

DEG_spleen_D22_qPCR_mod[DEG_spleen_D22_qPCR_mod$SYMBOL == target_qPCR_gene & grepl("Con_N", DEG_spleen_D22_qPCR_mod$contrast),]
sum_per_gene_tables_actin_2_5x[sum_per_gene_tables_actin_2_5x$SYMBOL == target_qPCR_gene,]
sum_per_gene_tables_actin_5x[sum_per_gene_tables_actin_5x$SYMBOL == target_qPCR_gene,]
sum_per_gene_tables_actin_10x[sum_per_gene_tables_actin_10x$SYMBOL == target_qPCR_gene,]

DEG_spleen_D22_qPCR_mod[DEG_spleen_D22_qPCR_mod$SYMBOL == target_qPCR_gene & grepl("Con_N", DEG_spleen_D22_qPCR_mod$contrast),]
sum_per_gene_tables_gapdh_2_5x[sum_per_gene_tables_gapdh_2_5x$SYMBOL == target_qPCR_gene,]
sum_per_gene_tables_gapdh_5x[sum_per_gene_tables_gapdh_5x$SYMBOL == target_qPCR_gene,]
sum_per_gene_tables_gapdh_10x[sum_per_gene_tables_gapdh_10x$SYMBOL == target_qPCR_gene,]

## filter the matched genes/condition in gene_condition_summary_ to get the final plot

passed_qPCR_genes <- c("HSPH1","DNAJA4","NEDD4L","ZBTB32","PCK1","CKMT2","ACSBG1","UBAP2L2","ASB1")
length(passed_qPCR_genes)

passed_qPCR_genes_actin_200_2_5x <- c("NEDD4L")
passed_qPCR_genes_actin_200_5x <- c("ACSBG1")
passed_qPCR_genes_actin_200_10x <- c("")

passed_qPCR_genes_actin_400_2_5x <- c("HSPH1","ZBTB32","PCK1","CKMT2","ASB1")
passed_qPCR_genes_actin_400_5x <- c("DNAJA4","UBAP2L2")

## run on 1/200 mode

passed_gene_condition_summary_actin_200_2_5x <- gene_condition_summary_actin_2_5x %>%
  filter(target_gene %in% passed_qPCR_genes_actin_200_2_5x)
passed_gene_condition_summary_actin_200_5x <- gene_condition_summary_actin_5x %>%
  filter(target_gene %in% passed_qPCR_genes_actin_200_5x)

## run on 1/400 mode

passed_gene_condition_summary_actin_400_2_5x <- gene_condition_summary_actin_2_5x %>%
  filter(target_gene %in% passed_qPCR_genes_actin_400_2_5x)
passed_gene_condition_summary_actin_400_5x <- gene_condition_summary_actin_5x %>%
  filter(target_gene %in% passed_qPCR_genes_actin_400_5x)

passed_gene_condition_summary_actin <- rbind(passed_gene_condition_summary_actin_200_2_5x,
                                             passed_gene_condition_summary_actin_200_5x,
                                             passed_gene_condition_summary_actin_400_2_5x,
                                             passed_gene_condition_summary_actin_400_5x)

## plot the final fig

## plot parameters

grey_map <- c(
  Con_N   = "black",
  TM_N    = "grey25",
  Con_AHS = "grey40",
  TM_AHS  = "grey55", 
  Con_CHS = "grey45",
  TM_CHS  = "grey65")

passed_gene_condition_summary_actin <- passed_gene_condition_summary_actin %>%
  mutate(target_gene = paste0("italic('", target_gene, "')"))

passed_qPCR_genes_levels <- c("italic('NEDD4L')",
                              "italic('ZBTB32')",
                              "italic('UBAP2L2')",
                              "italic('HSPH1')",
                              "italic('DNAJA4')",
                              "italic('CKMT2')",
                              "italic('PCK1')",
                              "italic('ASB1')",
                              "italic('ACSBG1')"
)

passed_qPCR_genes_levels <- c("italic('NEDD4L')",
                              "italic('ZBTB32')",
                              "italic('UBAP2L2')",
                              "italic('CKMT2')",
                              "italic('HSPH1')",
                              "italic('DNAJA4')",
                              "italic('PCK1')",
                              "italic('ASB1')"
                              # "italic('ACSBG1')"
)

passed_gene_condition_summary_actin$target_gene <- factor(passed_gene_condition_summary_actin$target_gene,
                                                          levels=passed_qPCR_genes_levels)
passed_gene_condition_summary_actin <- passed_gene_condition_summary_actin[order(passed_gene_condition_summary_actin$target_gene), ]

passed_qPCR_genes_levels <- list(
  set1 = c("italic('NEDD4L')", "italic('ZBTB32')","italic('UBAP2L2')"),
  set2 = c("italic('HSPH1')", "italic('DNAJA4')","italic('CKMT2')"),
  set3 = c("italic('PCK1')", "italic('ASB1')", "italic('ACSBG1')")
  )

## final plot

ggplot(subset(passed_gene_condition_summary_actin, 
              # target_gene != "ACSBG1"
              target_gene %in% passed_qPCR_genes_levels$set3
),
aes(x = condition_label, y = mean_fold_change,
    fill = condition_label, color = condition_label)) +
  geom_col(width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_fold_change - se_fold_change,
        ymax = mean_fold_change + se_fold_change),
    width = 0.2, color = "black") +
  geom_text(
    aes(y = annotation_y_pos, label = significance_symbol),
    vjust = 0, color = "black") +
  scale_fill_manual(values = grey_map) +
  scale_color_manual(values = grey_map) +
  facet_wrap(~ target_gene, scales = "free", ncol = 3, labeller = label_parsed) +
  labs(x = NULL, y = "Relative expression (2^-ΔΔCt)") +
  theme_classic() +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none", 
    strip.background = element_blank(),
    strip.text = element_text(face = "bold",size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1))

## compare actin and gapdh HK genes results

# no need to filter for passed; instead test for all genes after removing NULL
# passed_qPCR_genes_actin_2_5x <- c("HSPH1","DNAJA4","LOC107049986","GADD45G")
# passed_qPCR_genes_actin_5x <- c("COL15A1")

sum_per_gene_tables_actin_2_5x_mod <- sum_per_gene_tables_actin_2_5x %>% drop_na()
sum_per_gene_tables_gapdh_2_5x_mod <- sum_per_gene_tables_gapdh_2_5x %>% drop_na()
sum_per_gene_tables_actin_5x_mod <- sum_per_gene_tables_actin_5x %>% drop_na()
sum_per_gene_tables_gapdh_5x_mod <- sum_per_gene_tables_gapdh_5x %>% drop_na()
sum_per_gene_tables_actin_10x_mod <- sum_per_gene_tables_actin_10x %>% drop_na()
sum_per_gene_tables_gapdh_10x_mod <- sum_per_gene_tables_gapdh_10x %>% drop_na()

nrow(sum_per_gene_tables_actin_2_5x)
nrow(sum_per_gene_tables_gapdh_2_5x)
nrow(sum_per_gene_tables_actin_5x)
nrow(sum_per_gene_tables_gapdh_5x)
nrow(sum_per_gene_tables_actin_10x)
nrow(sum_per_gene_tables_gapdh_10x)

nrow(sum_per_gene_tables_actin_2_5x_mod)
nrow(sum_per_gene_tables_gapdh_2_5x_mod)
nrow(sum_per_gene_tables_actin_5x_mod)
nrow(sum_per_gene_tables_gapdh_5x_mod)
nrow(sum_per_gene_tables_actin_10x_mod)
nrow(sum_per_gene_tables_gapdh_10x_mod)

sum_per_gene_tables_actin_2_5x_mod$housekeeping <- "ACTB"
sum_per_gene_tables_gapdh_2_5x_mod$housekeeping <- "GAPDH"
sum_per_gene_tables_actin_5x_mod$housekeeping <- "ACTB"
sum_per_gene_tables_gapdh_5x_mod$housekeeping <- "GAPDH"
sum_per_gene_tables_actin_10x_mod$housekeeping <- "ACTB"
sum_per_gene_tables_gapdh_10x_mod$housekeeping <- "GAPDH"

actin_gapdh_comp_2_5x <- rbind(sum_per_gene_tables_actin_2_5x_mod, sum_per_gene_tables_gapdh_2_5x_mod)
actin_gapdh_comp_5x <- rbind(sum_per_gene_tables_actin_5x_mod, sum_per_gene_tables_gapdh_5x_mod)
actin_gapdh_comp_10x <- rbind(sum_per_gene_tables_actin_10x_mod, sum_per_gene_tables_gapdh_10x_mod)

hk_comp_2_5x <- actin_gapdh_comp_2_5x %>%
  select(SYMBOL, Name, fc, housekeeping) %>%
  tidyr::pivot_wider(names_from = housekeeping, values_from = fc)
hk_comp_5x <- actin_gapdh_comp_5x %>%
  select(SYMBOL, Name, fc, housekeeping) %>%
  tidyr::pivot_wider(names_from = housekeeping, values_from = fc)
hk_comp_10x <- actin_gapdh_comp_10x %>%
  select(SYMBOL, Name, fc, housekeeping) %>%
  tidyr::pivot_wider(names_from = housekeeping, values_from = fc)

hk_comp_2_5x$diff <- round((hk_comp_2_5x$ACTB - hk_comp_2_5x$GAPDH),7)
hk_comp_5x$diff <- round((hk_comp_5x$ACTB - hk_comp_5x$GAPDH),7)
hk_comp_10x$diff <- round((hk_comp_10x$ACTB - hk_comp_10x$GAPDH),7)

summary(hk_comp_2_5x$diff)
summary(hk_comp_5x$diff)
summary(hk_comp_10x$diff)

cor(hk_comp_2_5x$ACTB, hk_comp_2_5x$GAPDH, use = "complete.obs")
cor(hk_comp_5x$ACTB, hk_comp_5x$GAPDH, use = "complete.obs")
cor(hk_comp_10x$ACTB, hk_comp_10x$GAPDH, use = "complete.obs")

plot(log2(hk_comp_2_5x$ACTB),log2(hk_comp_2_5x$GAPDH), pch = 16)
abline(lm(log2(hk_comp_2_5x$GAPDH) ~ log2(hk_comp_2_5x$ACTB)), col = "firebrick", lwd = 2)
plot(log2(hk_comp_5x$ACTB),log2(hk_comp_5x$GAPDH), pch = 16)
abline(lm(log2(hk_comp_5x$GAPDH) ~ log2(hk_comp_5x$ACTB)), col = "firebrick", lwd = 2)
plot(log2(hk_comp_10x$ACTB),log2(hk_comp_10x$GAPDH), pch = 16)
abline(lm(log2(hk_comp_10x$GAPDH) ~ log2(hk_comp_10x$ACTB)), col = "firebrick", lwd = 2)

## quick and dirty model of X vs 2^-X

x <- seq(-50, 50, by = 0.5)
y <- 2^(-x)

model_expo <- data.frame(X = x, Y = y)

ggplot(model_expo, aes(x = X, y = Y)) +
  geom_line(size = 1) +
  labs(x = "X", y = expression(2^{-X})) +
  theme_classic()

## Compare with RNA-seq results - MultTissue_Devo   [in progress] ======================

## read DEG RNA-seq results

DEG_MultTissue_Devo_qPCR <- read_excel("C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_MultTissue_Devo\\TMomics_DEG_MultiTissueDevo_ctrl_pertissue_lfcShrink_apeglm_allanno_allcomps_final.xlsx", 
                                  sheet = "Sheet1")

unique(DEG_MultTissue_Devo_qPCR$contrast)
unique(DEG_MultTissue_Devo_qPCR$SYMBOL)
nrow(DEG_MultTissue_Devo_qPCR)

DEG_MultTissue_Devo_qPCR_mod <- DEG_MultTissue_Devo_qPCR %>%
  select(c("SYMBOL","log2FoldChange","direction","contrast")) %>%
  drop_na()
  # filter(grepl("liver", contrast, ignore.case = TRUE))

unique(DEG_MultTissue_Devo_qPCR_mod$SYMBOL)
nrow(DEG_MultTissue_Devo_qPCR_mod)
unique(DEG_MultTissue_Devo_qPCR_mod$contrast)

qPCR_genes <- c(
  "THRSP", "CYP7A1", "PIGR", "SCD",
  "TPM2-149", "COL9A1-199", "WNT6-86", "MATN4-111", "CHRNA9-130",
  "JCHAIN-84", "SPIN1L-119", "AVBD13-87", "LOC-146", "SULT-124"
)
length(unique(DEG_MultTissue_Devo_qPCR_mod$SYMBOL))
length(qPCR_genes)

## read qPCR processed results

base_dir <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\qPCR_MultTissue_Devo\\Liver\\"
read_per_gene_file <- function(path) {
  shs <- readxl::excel_sheets(path)
  dplyr::bind_rows(lapply(shs, function(s) {
    df <- readxl::read_excel(path, sheet = s)
    df$SYMBOL <- s
    df
  }))
}

## Liver
f_actin_2_5x <- file.path(base_dir, "Liver_Devo_CT1_CT2_cDNA_Final_SS_mod_ACTIN_2_5x_per_gene_tables.xlsx")
f_gapdh_2_5x <- file.path(base_dir, "Liver_Devo_CT1_CT2_cDNA_Final_SS_mod_GAPDH_2_5x_per_gene_tables.xlsx")

# ## Muscle
# f_actin_2_5x <- file.path(base_dir, "Muscle_Devo_CT1_CT2_cDNA_Final_SS_mod_ACTIN_2_5x_per_gene_tables.xlsx")
# f_gapdh_2_5x <- file.path(base_dir, "Muscle_Devo_CT1_CT2_cDNA_Final_SS_mod_GAPDH_2_5x_per_gene_tables.xlsx")

# ## Spleen
# f_actin_2_5x <- file.path(base_dir, "Spleen_Devo_CT1_CT2_cDNA_Final_SS_mod_ACTIN_2_5x_per_gene_tables.xlsx")
# f_gapdh_2_5x <- file.path(base_dir, "Spleen_Devo_CT1_CT2_cDNA_Final_SS_mod_GAPDH_2_5x_per_gene_tables.xlsx")

per_gene_tables_actin_2_5x <- read_per_gene_file(f_actin_2_5x)
per_gene_tables_gapdh_2_5x <- read_per_gene_file(f_gapdh_2_5x)

sum_per_gene_tables_actin_2_5x <- summarise(group_by(per_gene_tables_actin_2_5x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")
sum_per_gene_tables_gapdh_2_5x <- summarise(group_by(per_gene_tables_gapdh_2_5x,SYMBOL,Name), fc = mean(`2^-ΔΔCt`, na.rm = TRUE), .groups = "drop")

cond_levels <- c("D19Liver", "D7Liver", "D22Liver")
# cond_levels <- c("D19Muscle", "D7Muscle", "D22Muscle")
# cond_levels <- c("D19Spleen", "D7Spleen", "D22Spleen")

sum_per_gene_tables_actin_2_5x$Name <- factor(sum_per_gene_tables_actin_2_5x$Name,levels=cond_levels)
sum_per_gene_tables_gapdh_2_5x$Name <- factor(sum_per_gene_tables_gapdh_2_5x$Name,levels=cond_levels)

sum_per_gene_tables_actin_2_5x <- sum_per_gene_tables_actin_2_5x[order(sum_per_gene_tables_actin_2_5x$Name), ]
sum_per_gene_tables_gapdh_2_5x <- sum_per_gene_tables_gapdh_2_5x[order(sum_per_gene_tables_gapdh_2_5x$Name), ]

## go through each gene_condition_summary_ or sum_per_gene_tables_ to compare expression pattern with DEG_MultTissue_Devo_qPCR_mod

## (1) gene_condition_summary_actin_2_5x    OR    sum_per_gene_tables_actin_2_5x
## (2) gene_condition_summary_gapdh_2_5x    OR    sum_per_gene_tables_gapdh_2_5x

qPCR_genes <- c(
  "THRSP", "CYP7A1", "PIGR", "SCD",
  "TPM2-149", "COL9A1-199", "WNT6-86", "MATN4-111", "CHRNA9-130",
  "JCHAIN-84", "SPIN1L-119", "AVBD13-87", "LOC-146", "SULT-124"
)

target_qPCR_gene        <- c("SULT") 
target_qPCR_gene_primer <- c("SULT-124") 

DEG_MultTissue_Devo_qPCR_mod[DEG_MultTissue_Devo_qPCR_mod$SYMBOL == target_qPCR_gene & grepl("D19", DEG_MultTissue_Devo_qPCR_mod$contrast),]
sum_per_gene_tables_actin_2_5x[sum_per_gene_tables_actin_2_5x$SYMBOL == target_qPCR_gene_primer,]
# sum_per_gene_tables_gapdh_2_5x[sum_per_gene_tables_gapdh_2_5x$SYMBOL == target_qPCR_gene,]

qPCR_genes_mod <- c(
  # "THRSP", "CYP7A1", "PIGR", "SCD",
  "TPM2", "COL9A1", "WNT6", "MATN4", "CHRNA9",
  "JCHAIN", "SPIN1L", "AVBD1", "LOC107049986", "SULT"
)

DEG_MultTissue_Devo_qPCR_mod[DEG_MultTissue_Devo_qPCR_mod$SYMBOL %in% qPCR_genes_mod,]

# note that in spleen LOC107049986	is also XM_046902442; u can use to check
# DEG_MultTissue_Devo_qPCR %>% filter(grepl("XM_046902442", REFSEQ))

## filter the matched genes/condition in gene_condition_summary_ to get the final plot

## to automate later
## save gene_condition_summary_actin_2_5x and sum_per_gene_tables_actin_2_5x for each tissue

# rm(list = setdiff(ls(), c(
#   "sum_per_gene_tables_actin_Liver",
#   "sum_per_gene_tables_actin_Muscle",
#   "sum_per_gene_tables_actin_Spleen",
#   "gene_condition_summary_actin_Liver",
#   "gene_condition_summary_actin_Muscle",
#   "gene_condition_summary_actin_Spleen"
#   )))

sum_per_gene_tables_actin_Liver  <- sum_per_gene_tables_actin_2_5x
sum_per_gene_tables_actin_Muscle <- sum_per_gene_tables_actin_2_5x
sum_per_gene_tables_actin_Spleen <- sum_per_gene_tables_actin_2_5x

gene_condition_summary_actin_Liver  <- gene_condition_summary_actin_2_5x
gene_condition_summary_actin_Muscle <- gene_condition_summary_actin_2_5x
gene_condition_summary_actin_Spleen <- gene_condition_summary_actin_2_5x

DEG_MultTissue_Devo_qPCR_mod

passed_qPCR_genes <- c("CYP7A1","PIGR","SCD","COL9A1-199","MATN4-111","JCHAIN-84","AVBD13-87","SULT-124")
length(passed_qPCR_genes)

passed_qPCR_genes_actin_Liver <- c("CYP7A1","PIGR","SCD")
passed_qPCR_genes_actin_Muscle <- c("COL9A1-199","MATN4-111")
passed_qPCR_genes_actin_Spleen <- c("JCHAIN-84","AVBD13-87","SULT-124")

## run on three tissues

passed_gene_condition_summary_actin_Liver <- gene_condition_summary_actin_Liver %>%
  filter(target_gene %in% passed_qPCR_genes_actin_Liver)
passed_gene_condition_summary_actin_Muscle <- gene_condition_summary_actin_Muscle %>%
  filter(target_gene %in% passed_qPCR_genes_actin_Muscle)
passed_gene_condition_summary_actin_Spleen <- gene_condition_summary_actin_Spleen %>%
  filter(target_gene %in% passed_qPCR_genes_actin_Spleen)

passed_gene_condition_summary_actin_AllTissues <- rbind(
  passed_gene_condition_summary_actin_Liver,
  passed_gene_condition_summary_actin_Muscle,
  passed_gene_condition_summary_actin_Spleen)

## plot the final fig

## plot parameters

grey_map <- c(
  D19Liver    = "black",
  D19Muscle   = "black",
  D19Spleen   = "black",
  D7Liver     = "grey25",
  D7Muscle    = "grey25",
  D7Spleen    = "grey25",
  D22Liver   = "grey65",
  D22Muscle  = "grey65",
  D22Spleen  = "grey65"
  )

grey_map <- c(
  ED19   = "black",
  D7    = "grey25",
  D22   = "grey65"
)

passed_gene_condition_summary_actin_AllTissues <- passed_gene_condition_summary_actin_AllTissues %>%
  mutate(target_gene = paste0("italic('", target_gene, "')"))

passed_qPCR_genes_levels <- c("italic('CYP7A1')",
                              "italic('PIGR')",
                              "italic('SCD')",
                              "italic('COL9A1-199')",
                              "italic('MATN4-111')",
                              "italic('JCHAIN-84')",
                              "italic('AVBD13-87')",
                              "italic('SULT-124')"
                              )

passed_gene_condition_summary_actin_AllTissues$target_gene <- factor(
  passed_gene_condition_summary_actin_AllTissues$target_gene,
  levels=passed_qPCR_genes_levels)
passed_gene_condition_summary_actin_AllTissues <- passed_gene_condition_summary_actin_AllTissues[order(
  passed_gene_condition_summary_actin_AllTissues$target_gene), ]

passed_qPCR_genes_levels <- list(
  set1 = c("italic('CYP7A1')", "italic('PIGR')","italic('SCD')"),
  set2 = c("italic('COL9A1-199')", "italic('MATN4-111')"),
  set3 = c("italic('JCHAIN-84')", "italic('AVBD13-87')", "italic('SULT-124')")
)

passed_gene_condition_summary_actin_AllTissues_mod <- passed_gene_condition_summary_actin_AllTissues %>%
  mutate(condition_label = case_when(
    grepl("D19", condition_label) ~ "ED19",
    grepl("D7", condition_label)  ~ "D7",
    grepl("D22", condition_label) ~ "D22",
    TRUE ~ condition_label
  ))  %>%
  mutate(condition_label = factor(condition_label, levels = c("ED19", "D7", "D22")))

## final plot

ggplot(subset(passed_gene_condition_summary_actin_AllTissues_mod, 
              # target_gene != "ACSBG1"
              # target_gene %in% passed_qPCR_genes_levels$set3
),
aes(x = condition_label, y = mean_fold_change,
    fill = condition_label, color = condition_label)) +
  geom_col(width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_fold_change - se_fold_change,
        ymax = mean_fold_change + se_fold_change),
    width = 0.2, color = "black") +
  geom_text(
    aes(y = annotation_y_pos, label = significance_symbol),
    vjust = 0, color = "black") +
  scale_fill_manual(values = grey_map) +
  scale_color_manual(values = grey_map) +
  facet_wrap(~ target_gene, scales = "free", ncol = 4, labeller = label_parsed) +
  labs(x = NULL, y = "Relative expression (2^-ΔΔCt)") +
  theme_classic() +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none", 
    strip.background = element_blank(),
    strip.text = element_text(face = "bold",size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 1))

## compare actin and gapdh HK genes results

# no need to filter for passed; instead test for all genes after removing NULL

sum_per_gene_tables_actin_2_5x_mod <- sum_per_gene_tables_actin_2_5x %>% drop_na()
sum_per_gene_tables_gapdh_2_5x_mod <- sum_per_gene_tables_gapdh_2_5x %>% drop_na()

nrow(sum_per_gene_tables_actin_2_5x)
nrow(sum_per_gene_tables_gapdh_2_5x)

nrow(sum_per_gene_tables_actin_2_5x_mod)
nrow(sum_per_gene_tables_gapdh_2_5x_mod)

sum_per_gene_tables_actin_2_5x_mod$housekeeping <- "ACTB"
sum_per_gene_tables_gapdh_2_5x_mod$housekeeping <- "GAPDH"

actin_gapdh_comp_2_5x <- rbind(sum_per_gene_tables_actin_2_5x_mod, sum_per_gene_tables_gapdh_2_5x_mod)

hk_comp_2_5x <- actin_gapdh_comp_2_5x %>%
  select(SYMBOL, Name, fc, housekeeping) %>%
  tidyr::pivot_wider(names_from = housekeeping, values_from = fc)

hk_comp_2_5x$diff <- round((hk_comp_2_5x$ACTB - hk_comp_2_5x$GAPDH),7)

summary(hk_comp_2_5x$diff)

cor(hk_comp_2_5x$ACTB, hk_comp_2_5x$GAPDH, use = "complete.obs")

plot(log2(hk_comp_2_5x$ACTB),log2(hk_comp_2_5x$GAPDH), pch = 16)
abline(lm(log2(hk_comp_2_5x$GAPDH) ~ log2(hk_comp_2_5x$ACTB)), col = "firebrick", lwd = 2)

## Combine RT-qPCR (2^-ΔΔCt) and RNA-seq (log2FoldChange) results
## into a single comparison plot
## contrast format: "D22 vs D19 | liver"  /  "D7 vs D19 | liver"
## D19 is the baseline -> placeholder log2FC = 0

## gene symbol mapping between qPCR primer names and RNA-seq SYMBOLs

qPCR_to_RNAseq_map <- c(
  "CYP7A1"      = "CYP7A1",
  "PIGR"        = "PIGR",
  "SCD"         = "SCD",
  "COL9A1-199"  = "COL9A1",
  "MATN4-111"   = "MATN4",
  "JCHAIN-84"   = "JCHAIN",
  "AVBD13-87"   = "AVBD1",
  "SULT-124"    = "SULT"
)

## tissue mapping per gene (from your passed_qPCR_genes_actin_* lists)

qPCR_gene_tissue_map <- c(
  "CYP7A1"      = "Liver",
  "PIGR"        = "Liver",
  "SCD"         = "Liver",
  "COL9A1-199"  = "Muscle",
  "MATN4-111"   = "Muscle",
  "JCHAIN-84"   = "Spleen",
  "AVBD13-87"   = "Spleen",
  "SULT-124"    = "Spleen"
)

## prep qPCR side 

qPCR_plot_df <- passed_gene_condition_summary_actin_AllTissues_mod %>%
  mutate(target_gene_clean = gsub("italic\\('|'\\)", "", target_gene)) %>%
  mutate(RNAseq_SYMBOL = qPCR_to_RNAseq_map[target_gene_clean],
         tissue        = qPCR_gene_tissue_map[target_gene_clean]) %>%
  select(SYMBOL = RNAseq_SYMBOL, tissue, condition_label,
         value = mean_fold_change, se = se_fold_change,
         significance_symbol) %>%
  mutate(value_log2 = log2(value)) %>%
  mutate(method = "RT-qPCR")

## prep RNA-seq side
## split contrast into comparison + tissue, keep only *vs D19* comparisons

RNAseq_plot_df <- DEG_MultTissue_Devo_qPCR_mod %>%
  filter(SYMBOL %in% qPCR_plot_df$SYMBOL) %>%
  tidyr::separate(contrast, into = c("comparison", "tissue_raw"),
                  sep = "\\|", remove = FALSE) %>%
  mutate(comparison = trimws(comparison),
         tissue_raw = trimws(tissue_raw)) %>%
  filter(comparison %in% c("D22 vs D19", "D7 vs D19")) %>%
  mutate(tissue = case_when(
    grepl("liver",  tissue_raw, ignore.case = TRUE) ~ "Liver",
    grepl("muscle", tissue_raw, ignore.case = TRUE) ~ "Muscle",
    grepl("spleen", tissue_raw, ignore.case = TRUE) ~ "Spleen",
    TRUE ~ NA_character_
  )) %>%
  mutate(condition_label = case_when(
    comparison == "D22 vs D19" ~ "D22",
    comparison == "D7 vs D19"  ~ "D7",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(tissue), !is.na(condition_label)) %>%
  inner_join(
    qPCR_plot_df %>% distinct(SYMBOL, tissue),
    by = c("SYMBOL", "tissue")
  ) %>%
  select(SYMBOL, tissue, condition_label, value_log2 = log2FoldChange) %>%
  mutate(value = 2^value_log2,
         se = NA_real_,
         significance_symbol = NA_character_,
         method = "RNA-seq")

## add D19 baseline placeholder rows for RNA-seq (log2FC = 0 by definition)

RNAseq_D19_placeholder <- RNAseq_plot_df %>%
  distinct(SYMBOL, tissue) %>%
  mutate(condition_label = "ED19",
         value_log2 = 0,
         value = 1,
         se = NA_real_,
         significance_symbol = NA_character_,
         method = "RNA-seq")

RNAseq_plot_df <- bind_rows(RNAseq_plot_df, RNAseq_D19_placeholder)

## combine

combined_plot_df <- bind_rows(
  qPCR_plot_df  %>% select(SYMBOL, tissue, condition_label, value, value_log2, se, significance_symbol, method),
  RNAseq_plot_df %>% select(SYMBOL, tissue, condition_label, value, value_log2, se, significance_symbol, method)
)

combined_plot_df$condition_label <- factor(
  combined_plot_df$condition_label, levels = c("ED19", "D7", "D22"))

combined_plot_df$method <- factor(
  combined_plot_df$method, levels = c("RT-qPCR", "RNA-seq"))

combined_plot_df <- combined_plot_df %>%
  mutate(SYMBOL_label = paste0("italic('", SYMBOL, "')"))

symbol_levels <- paste0("italic('", c(
  "CYP7A1","PIGR","SCD","COL9A1","MATN4","JCHAIN","AVBD1","SULT"), "')")

combined_plot_df$SYMBOL_label <- factor(
  combined_plot_df$SYMBOL_label, levels = symbol_levels)

## plot parameters

method_map <- c("RT-qPCR" = "grey35", "RNA-seq" = "steelblue4")

## final combined plot

ggplot(combined_plot_df,
       aes(x = condition_label, y = value_log2,
           fill = method, group = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black") +
  geom_errorbar(
    data = subset(combined_plot_df, method == "RT-qPCR"),
    aes(ymin = log2(value - se), ymax = log2(value + se)),
    position = position_dodge(width = 0.7), width = 0.2, color = "black") +
  geom_text(
    data = subset(combined_plot_df, method == "RT-qPCR"),
    aes(y = log2(value) + 0.3, label = significance_symbol),
    position = position_dodge(width = 0.7), vjust = 0, color = "black") +
  scale_fill_manual(values = method_map, name = NULL) +
  facet_wrap(~ SYMBOL_label, scales = "free", ncol = 3, labeller = label_parsed) +
  labs(x = NULL, y = expression(log[2]~"(fold change vs D19)")) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5))

## correlation between RT-qPCR and RNA-seq log2FC per gene/condition
## (D19 excluded since it's a fixed placeholder with no variance)

corr_df <- combined_plot_df %>%
  filter(condition_label != "ED19") %>%
  select(SYMBOL, tissue, condition_label, method, value_log2) %>%
  tidyr::pivot_wider(names_from = method, values_from = value_log2)

cor(corr_df$`RT-qPCR`, corr_df$`RNA-seq`, use = "complete.obs")

## Pearson correlation

cor_test <- cor.test(
  corr_df$`RT-qPCR`,
  corr_df$`RNA-seq`,
  use = "complete.obs",
  method = "pearson")

r_value <- round(cor_test$estimate, 4)
p_value <- round(cor_test$p.value, 4)

ggplot(
  corr_df,
  aes(x = `RNA-seq`, y = `RT-qPCR`, label = SYMBOL)) +
  geom_point(size = 2.5, color = "grey20") +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "firebrick",
    linewidth = 0.7) +
  ggrepel::geom_text_repel(size = 4, fontface = "italic") +
  annotate("text", x = -Inf, y = Inf, label = paste0("R = ", r_value, "\n", "p = ", p_value),
           hjust = -0.1, vjust = 1.5, size = 4.5)+
  labs(
    x = expression("RNA-seq log"[2]*"FC"),
    y = expression("RT-qPCR log"[2]*"FC")) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(fill = NA, linewidth = 0.8),
    axis.text = element_text(color = "black"))

## Compare with RNA-seq results - Atro_Ex           [in progress] --------------------


