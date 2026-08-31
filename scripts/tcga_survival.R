# ============================================================
# tcga_survival.R
# Associate HPV pathway gene sets (from spatial GO/KEGG enrichment)
# with TCGA-HNSC bulk RNA-seq survival outcomes.
#
# Inputs:
#   go_results     — output of run_go()  from st_hpv_pathway_analysis.R
#   cache_folder   — directory for .rds checkpoints (shared with main pipeline)
#   picture_folder — directory for PNG/CSV outputs
#   clinical_csv   — (optional) path to a pre-processed TCGA clinical CSV
#                    (survival_data_HNSC.csv). When supplied, this is used
#                    instead of downloading clinical data via TCGAbiolinks.
#                    Expected columns: cases.submitter_id, overall_survival,
#                    deceased (TRUE/FALSE), Subtype (HNSC_HPV+ / HNSC_HPV-).
#
# Workflow:
#   1. Download TCGA-HNSC RNA-seq counts + clinical data  (cached)
#      OR load pre-processed clinical CSV
#   2. Extract gene sets from spatial GO enrichment results
#   3. Score gene sets in bulk RNA-seq with GSVA
#   4. Kaplan-Meier curves — high vs low score, optionally split by HPV status
#   5. Multivariate Cox proportional hazards regression
# ============================================================

library(TCGAbiolinks)
library(SummarizedExperiment)
library(GSVA)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(tibble)
library(GEOquery)      # NEW — for fetching GEO cohorts


# ── Config ─────────────────────────────────────────────────────────────────────

# Minimum number of genes required in a set to run GSVA / Cox
MIN_GENESET_SIZE <- 5

# Score dichotomization: split each pathway score at this quantile
# (0.5 = median split; 0.67 = top vs bottom tertile, etc.)
SCORE_SPLIT_QUANTILE <- 0.5

# FDR cutoff for Cox results table
COX_FDR_CUTOFF <- 0.05

# Number of top GO terms to pull from each HPV direction
N_GO_TERMS <- 5

VALIDATION_COHORTS <- list(
  GSE65858 = list(
    accession    = "GSE65858",
    surv_time    = "overall_survival_month:ch1",   # months → converted below
    surv_event   = "vital_status:ch1",              # "dead" / "alive"
    event_level  = "dead",                          # which level = event
    time_unit    = "months",                        # "months" or "days"
    hpv_col      = "hpv_status:ch1"
  ),
  GSE41613 = list(
    accession    = "GSE41613",
    surv_time    = "disease-specific survival (days):ch1",
    surv_event   = "vital status:ch1",
    event_level  = "dead",
    time_unit    = "days",
    hpv_col      = NULL                            # no HPV in this cohort
  ),
  GSE42743 = list(
    accession    = "GSE42743",
    surv_time    = "overall survival (months):ch1",
    surv_event   = "vital status:ch1",
    event_level  = "deceased",
    time_unit    = "months",
    hpv_col      = NULL                            # all HPV+ by design
  )
)

# ── Step 1: Download / cache TCGA-HNSC data ───────────────────────────────────

#' Download TCGA-HNSC RNA-seq (STAR counts) and clinical data, cache as .rds.
#'
#' @param cache_path  Directory to store cached files.
#' @return list(expr = SummarizedExperiment, clinical = data.frame)
get_tcga_hnsc <- function(cache_path) {
  tcga_cache <- file.path(cache_path, "tcga_hnsc.rds")

  if (file.exists(tcga_cache)) {
    cat("Loading TCGA-HNSC from cache...\n")
    return(readRDS(tcga_cache))
  }

  cat("Querying TCGA-HNSC RNA-seq data via TCGAbiolinks...\n")
  query_rna <- GDCquery(
    project           = "TCGA-HNSC",
    data.category     = "Transcriptome Profiling",
    data.type         = "Gene Expression Quantification",
    workflow.type     = "STAR - Counts"
  )

  cat("Downloading files (this may take several minutes)...\n")
  GDCdownload(query_rna, method = "api", files.per.chunk = 10)

  cat("Preparing SummarizedExperiment...\n")
  se <- GDCprepare(query_rna)

  cat("Downloading clinical data...\n")
  clinical <- GDCquery_clinic("TCGA-HNSC", type = "clinical")

  cat("Downloading TCGA-HNSC subtype data (HPV status)...\n")
  subtype <- TCGAquery_subtype(tumor = "HNSC")

  result <- list(expr = se, clinical = clinical, subtype = subtype)
  dir.create(cache_path, recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, tcga_cache)
  cat("TCGA-HNSC data cached to:", tcga_cache, "\n")
  return(result)
}


# ── Step 2: Tidy clinical data ─────────────────────────────────────────────────

#' Extract and clean survival + HPV status from TCGA clinical data.
#'
#' @param clinical  data.frame from GDCquery_clinic()
#' @param subtype   data.frame from TCGAquery_subtype("HNSC"), or NULL.
#'                  Used to supply HPV status when it is absent from clinical.
#' @return data.frame with columns: submitter_id, os_time, os_event, hpv_status
tidy_clinical <- function(clinical, subtype = NULL) {
  cat("Tidying clinical data...\n")

  clin <- clinical |>
    select(
      submitter_id,
      days_to_death,
      days_to_last_follow_up,
      vital_status,
      any_of(c("hpv_status", "hpv_calls", "paper_HPV.status"))
    ) |>
    mutate(
      os_time = case_when(
        !is.na(days_to_death)          ~ as.numeric(days_to_death),
        !is.na(days_to_last_follow_up) ~ as.numeric(days_to_last_follow_up),
        TRUE                           ~ NA_real_
      ),
      os_event = ifelse(tolower(vital_status) == "dead", 1, 0)
    ) |>
    filter(!is.na(os_time), os_time > 0)

  # Normalise any inline HPV column to "hpv_status"
  hpv_col <- intersect(c("hpv_status", "hpv_calls", "paper_HPV.status"), colnames(clin))[1]
  if (!is.na(hpv_col)) {
    clin <- clin |> rename(hpv_status = all_of(hpv_col))
  } else if (!is.null(subtype) && "paper_HPV.status" %in% colnames(subtype)) {
    # Fall back to subtype table — join on 12-char patient barcode
    cat("HPV status not in clinical data — pulling from subtype table.\n")
    hpv_sub <- subtype |>
      select(patient, hpv_status = paper_HPV.status) |>
      distinct()
    clin <- left_join(clin, hpv_sub, by = c("submitter_id" = "patient"))
  } else {
    cat("Warning: no HPV status found — HPV-stratified plots will be skipped.\n")
    clin$hpv_status <- NA_character_
  }

  # Standardize to "positive" / "negative" / NA
  clin$hpv_status <- tolower(trimws(clin$hpv_status))
  clin$hpv_status <- case_when(
    grepl("pos|yes|\\+", clin$hpv_status) ~ "positive",
    grepl("neg|no|\\-",  clin$hpv_status) ~ "negative",
    TRUE                                   ~ NA_character_
  )

  cat("Clinical data rows retained:", nrow(clin), "\n")
  cat("HPV status breakdown:\n")
  print(table(clin$hpv_status, useNA = "always"))
  return(clin)
}


# ── Step 2b: Load pre-processed clinical CSV (preferred over TCGAbiolinks) ────

#' Load survival + HPV status from a pre-processed TCGA CSV.
#'
#' Produces the same output format as tidy_clinical() so the rest of the
#' pipeline is unaffected.
#'
#' Expected CSV columns:
#'   cases.submitter_id  — 12-char TCGA patient barcode
#'   overall_survival    — survival time in days (pre-computed)
#'   deceased            — TRUE/FALSE event indicator
#'   Subtype             — "HNSC_HPV+" or "HNSC_HPV-"
#'   diagnoses.days_to_last_follow_up — used as fallback when overall_survival is NA
#'
#' @param csv_path  Full path to the CSV file.
#' @return data.frame with columns: submitter_id, os_time, os_event, hpv_status
load_external_clinical <- function(csv_path) {
  cat("Loading pre-processed clinical data from:", csv_path, "\n")
  raw <- read.csv(csv_path, stringsAsFactors = FALSE)

  # ── Validate required columns ──────────────────────────────────────────────
  required <- c("cases.submitter_id", "overall_survival", "deceased", "Subtype")
  missing_cols <- setdiff(required, colnames(raw))
  if (length(missing_cols) > 0) {
    stop("CSV is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  clin <- raw |>
    dplyr::mutate(
      submitter_id = cases.submitter_id,

      # overall_survival is pre-computed; fall back to days_to_last_follow_up
      # for the rare case where days_to_death was stored as '--' in the source
      os_time = dplyr::case_when(
        !is.na(overall_survival) & overall_survival > 0 ~ as.numeric(overall_survival),
        "diagnoses.days_to_last_follow_up" %in% names(raw) ~
          suppressWarnings(as.numeric(diagnoses.days_to_last_follow_up)),
        TRUE ~ NA_real_
      ),

      os_event = as.integer(deceased == TRUE | deceased == "TRUE"),

      # Standardize "HNSC_HPV+" / "HNSC_HPV-" → "positive" / "negative"
      hpv_status = dplyr::case_when(
        grepl("HPV\\+", Subtype) ~ "positive",
        grepl("HPV-",  Subtype) ~ "negative",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(os_time), os_time > 0) |>
    dplyr::select(submitter_id, os_time, os_event, hpv_status)

  cat("Clinical data rows retained:", nrow(clin), "\n")
  cat("HPV status breakdown:\n")
  print(table(clin$hpv_status, useNA = "always"))
  return(clin)
}


# ── Step 3: Build normalized expression matrix ────────────────────────────────

#' Extract log-normalized (log2 TPM+1) expression matrix from SummarizedExperiment.
#' GSVA works on genes x samples matrices.
#'
#' @param se  SummarizedExperiment from GDCprepare()
#' @return matrix of log2(TPM+1), rows = gene symbols, cols = TCGA barcodes (short)
build_expr_matrix <- function(se) {
  cat("Building expression matrix...\n")

  # Use TPM if available, otherwise fall back to unstranded counts
  if ("tpm_unstrand" %in% assayNames(se)) {
    mat <- assay(se, "tpm_unstrand")
    mat <- log2(mat + 1)
    cat("Using log2(TPM+1)\n")
  } else {
    mat <- assay(se, "unstranded")
    # Crude CPM normalization if TPM not available
    mat <- sweep(mat, 2, colSums(mat), "/") * 1e6
    mat <- log2(mat + 1)
    cat("TPM not found — using log2(CPM+1)\n")
  }

  # Use gene symbols as row names (fall back to gene_id if no symbol)
  gene_info  <- as.data.frame(rowData(se))
  symbol_col <- intersect(c("gene_name", "external_gene_name"), colnames(gene_info))[1]
  if (!is.na(symbol_col)) {
    symbols    <- gene_info[[symbol_col]]
    # Drop rows with missing/duplicated symbols
    keep       <- !is.na(symbols) & symbols != "" & !duplicated(symbols)
    mat        <- mat[keep, , drop = FALSE]
    rownames(mat) <- symbols[keep]
  }

  # Shorten TCGA barcodes to 12-char patient ID (TCGA-XX-XXXX)
  colnames(mat) <- substr(colnames(mat), 1, 12)
  # Keep one column per patient (primary tumors; drop duplicates)
  mat <- mat[, !duplicated(colnames(mat)), drop = FALSE]

  cat("Expression matrix:", nrow(mat), "genes x", ncol(mat), "samples\n")
  return(mat)
}


# ── Step 4: Extract gene sets from spatial GO results ─────────────────────────

#' Pull top N GO gene sets from run_go() output.
#'
#' @param go_results  list(HPV_pos = enrichResult, HPV_neg = enrichResult)
#' @param n_terms     number of top terms to extract per direction
#' @return named list of character vectors (gene symbols)
extract_go_genesets <- function(go_results, n_terms = N_GO_TERMS) {
  extract_sets <- function(ego, prefix) {
    if (is.null(ego) || nrow(ego@result) == 0) return(list())
    top <- head(ego@result[ego@result$p.adjust < 0.05, ], n_terms)
    if (nrow(top) == 0) return(list())
    sets <- lapply(top$geneID, function(x) unlist(strsplit(x, "/")))
    # Sanitize names for column/file naming
    names(sets) <- paste0(prefix, "_", gsub("[^A-Za-z0-9]", "_", substr(top$Description, 1, 40)))
    return(sets)
  }

  pos_sets <- extract_sets(go_results$HPV_pos, "HPVpos")
  neg_sets <- extract_sets(go_results$HPV_neg, "HPVneg")

  all_sets <- c(pos_sets, neg_sets)
  # Drop sets too small for GSVA
  all_sets <- all_sets[sapply(all_sets, length) >= MIN_GENESET_SIZE]
  cat(length(all_sets), "gene sets extracted from GO results\n")
  return(all_sets)
}


# ── Step 5: GSVA pathway scoring ──────────────────────────────────────────────

#' Score gene sets in TCGA bulk expression using GSVA.
#'
#' @param expr_mat   log-normalized matrix (genes x samples)
#' @param gene_sets  named list of gene symbol vectors
#' @return data.frame: rows = samples, cols = pathway scores + submitter_id
run_gsva <- function(expr_mat, gene_sets) {
  cat("Running GSVA on", length(gene_sets), "gene sets x", ncol(expr_mat), "samples...\n")

  # Filter sets to genes present in expression matrix
  gene_sets_filt <- lapply(gene_sets, function(g) intersect(g, rownames(expr_mat)))
  gene_sets_filt <- gene_sets_filt[sapply(gene_sets_filt, length) >= MIN_GENESET_SIZE]

  if (length(gene_sets_filt) == 0) stop("No gene sets had enough overlapping genes with the expression matrix.")

  cat("Gene sets after filtering for matrix overlap:", length(gene_sets_filt), "\n")
  for (nm in names(gene_sets_filt)) {
    cat("  ", nm, "—", length(gene_sets_filt[[nm]]), "genes\n")
  }

  gsva_param <- gsvaParam(expr_mat, gene_sets_filt, minSize = MIN_GENESET_SIZE)
  gsva_scores <- gsva(gsva_param, verbose = FALSE)

  scores_df <- as.data.frame(t(gsva_scores))
  scores_df$submitter_id <- rownames(scores_df)
  return(scores_df)
}


# ── Step 6: Kaplan-Meier survival curves ──────────────────────────────────────

#' Run KM survival analysis for each pathway score.
#' Dichotomizes scores at SCORE_SPLIT_QUANTILE and plots KM curves.
#'
#' @param scores_df    output of run_gsva()
#' @param clin         output of tidy_clinical()
#' @param out_path     directory for PNG output
#' @param stratify_hpv logical — also produce HPV-stratified curves
run_km_survival <- function(scores_df, clin, out_path, stratify_hpv = TRUE) {
  cat("Running Kaplan-Meier survival analysis...\n")

  # Join scores with clinical data
  combined <- inner_join(clin, scores_df, by = "submitter_id")
  cat("Samples with matched survival + scores:", nrow(combined), "\n")

  pathway_cols <- setdiff(colnames(scores_df), "submitter_id")

  for (pw in pathway_cols) {
    cat("  KM for:", pw, "\n")
    threshold <- quantile(combined[[pw]], SCORE_SPLIT_QUANTILE, na.rm = TRUE)
    combined$score_group <- ifelse(combined[[pw]] >= threshold, "High", "Low")

    # ── Overall KM ────────────────────────────────────────────────────────────
    fit <- survfit(Surv(os_time, os_event) ~ score_group, data = combined)
    p_km <- ggsurvplot(
      fit,
      data          = combined,
      pval          = TRUE,
      conf.int      = TRUE,
      risk.table    = TRUE,
      palette       = c("#E64B35", "#4DBBD5"),
      title         = paste0("Overall Survival — ", gsub("_", " ", pw)),
      xlab          = "Time (days)",
      legend.labs   = c("High score", "Low score"),
      ggtheme       = theme_classic()
    )
    # ggsurvplot returns a list, not a ggplot — must use png/dev.off to save correctly
    png(file.path(out_path, paste0("km_", pw, ".png")),
        width = 8, height = 7, units = "in", res = 150)
    print(p_km)
    dev.off()

    # ── HPV-stratified KM ─────────────────────────────────────────────────────
    if (stratify_hpv && any(!is.na(combined$hpv_status))) {
      for (hpv in c("positive", "negative")) {
        sub <- combined[!is.na(combined$hpv_status) & combined$hpv_status == hpv, ]
        if (nrow(sub) < 10 || length(unique(sub$score_group)) < 2) next

        fit_hpv <- survfit(Surv(os_time, os_event) ~ score_group, data = sub)
        p_hpv <- ggsurvplot(
          fit_hpv,
          data        = sub,
          pval        = TRUE,
          conf.int    = TRUE,
          risk.table  = TRUE,
          palette     = c("#E64B35", "#4DBBD5"),
          title       = paste0("OS — ", gsub("_", " ", pw), " (HPV ", hpv, ")"),
          xlab        = "Time (days)",
          legend.labs = c("High score", "Low score"),
          ggtheme     = theme_classic()
        )
        png(file.path(out_path, paste0("km_", pw, "_HPV_", hpv, ".png")),
            width = 8, height = 7, units = "in", res = 150)
        print(p_hpv)
        dev.off()
      }
    }
  }
  cat("KM plots saved to:", out_path, "\n")
}


# ── Step 7: Cox proportional hazards regression ───────────────────────────────

#' Fit a Cox model for each pathway score (univariate) and one joint multivariate model.
#' Adjusts for HPV status and age where available.
#'
#' @param scores_df  output of run_gsva()
#' @param clin       output of tidy_clinical()
#' @param out_path   directory for CSV / PNG output
#' @return data.frame of univariate Cox results
run_cox_survival <- function(scores_df, clin, out_path) {
  cat("Running Cox proportional hazards regression...\n")

  combined <- inner_join(clin, scores_df, by = "submitter_id")
  pathway_cols <- setdiff(colnames(scores_df), "submitter_id")

  # ── Univariate Cox: one model per pathway ─────────────────────────────────
  cox_results <- lapply(pathway_cols, function(pw) {
    formula <- as.formula(paste0("Surv(os_time, os_event) ~ `", pw, "`"))
    fit <- tryCatch(coxph(formula, data = combined), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    s <- summary(fit)
    data.frame(
      pathway  = pw,
      HR       = exp(coef(fit)),
      CI_lower = s$conf.int[, "lower .95"],
      CI_upper = s$conf.int[, "upper .95"],
      p_value  = s$coefficients[, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  })
  cox_results <- bind_rows(cox_results)
  cox_results$p_adj <- p.adjust(cox_results$p_value, method = "BH")
  cox_results <- cox_results[order(cox_results$p_adj), ]

  write.csv(cox_results, file.path(out_path, "cox_univariate_results.csv"), row.names = FALSE)
  cat("Univariate Cox results:\n")
  print(cox_results)

  # ── Forest plot of univariate results ────────────────────────────────────
  forest_df <- cox_results |>
    mutate(
      label     = gsub("_", " ", pathway),
      sig       = ifelse(p_adj < COX_FDR_CUTOFF, "FDR < 0.05", "NS")
    )

  p_forest <- ggplot(forest_df, aes(x = HR, y = reorder(label, -HR),
                                     xmin = CI_lower, xmax = CI_upper,
                                     color = sig)) +
    geom_pointrange(size = 0.6) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("FDR < 0.05" = "#E64B35", "NS" = "grey60")) +
    labs(
      title  = "Univariate Cox — Pathway Score vs Overall Survival (TCGA-HNSC)",
      x      = "Hazard Ratio (95% CI)",
      y      = NULL,
      color  = NULL
    ) +
    theme_classic() +
    theme(legend.position = "bottom")

  ggsave(file.path(out_path, "cox_forest_plot.png"),
         p_forest, width = 9, height = max(4, nrow(forest_df) * 0.5 + 2), dpi = 150)

  # ── Multivariate Cox: all pathways + HPV + age ────────────────────────────
  covariate_cols <- pathway_cols
  has_hpv <- any(!is.na(combined$hpv_status))
  has_age <- "age_at_diagnosis" %in% colnames(combined)

  if (has_hpv)  covariate_cols <- c(covariate_cols, "hpv_status")
  if (has_age)  covariate_cols <- c(covariate_cols, "age_at_diagnosis")

  multi_formula <- as.formula(
    paste0("Surv(os_time, os_event) ~ ",
           paste0("`", covariate_cols, "`", collapse = " + "))
  )
  multi_fit <- tryCatch(coxph(multi_formula, data = combined), error = function(e) {
    cat("Warning: multivariate Cox failed —", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(multi_fit)) {
    multi_summary <- as.data.frame(summary(multi_fit)$coefficients)
    multi_summary$variable <- rownames(multi_summary)
    write.csv(multi_summary, file.path(out_path, "cox_multivariate_results.csv"), row.names = FALSE)
    cat("Multivariate Cox results saved.\n")
  }

  return(invisible(cox_results))
}

# ── Load a GEO validation cohort ─────────────────────────────────────────────

#' Download a GEO cohort's expression matrix and clinical survival data.
#'
#' Returns the same shape as TCGA outputs so downstream functions are reused
#' without modification:
#'   expr_mat  — genes x samples log-normalized matrix
#'   clin      — data.frame(submitter_id, os_time [days], os_event, hpv_status)
#'
#' @param cohort_cfg  one entry from VALIDATION_COHORTS list
#' @param cache_folder directory for .rds checkpoints
#' @return list(expr_mat = matrix, clin = data.frame)
load_geo_cohort <- function(cohort_cfg, cache_folder) {
  acc        <- cohort_cfg$accession
  cache_path <- file.path(cache_folder, paste0(acc, ".rds"))
  
  if (file.exists(cache_path)) {
    cat("Loading", acc, "from cache...\n")
    return(readRDS(cache_path))
  }
  
  cat("Downloading", acc, "from GEO...\n")
  gse <- getGEO(acc, GSEMatrix = TRUE, destdir = cache_folder)
  # getGEO returns a list when there are multiple platforms
  if (is.list(gse)) gse <- gse[[1]]
  
  # ── Expression matrix ──────────────────────────────────────────────────────
  expr_mat <- exprs(gse)
  
  # Map probe IDs → gene symbols via feature data
  fdata      <- fData(gse)
  symbol_col <- intersect(c("Gene Symbol", "GENE_SYMBOL", "gene_assignment",
                            "Symbol", "gene_symbol"), colnames(fdata))[1]
  
  if (!is.na(symbol_col)) {
    symbols <- fdata[[symbol_col]]
    # Some platforms store "GENE // description // ..." — take first token
    symbols <- sapply(strsplit(as.character(symbols), " // "), `[`, 1)
    symbols <- trimws(symbols)
    keep    <- !is.na(symbols) & symbols != "" & symbols != "---"
    
    expr_mat <- expr_mat[keep, , drop = FALSE]
    symbols  <- symbols[keep]
    
    # Collapse probes to gene level by keeping the probe with highest mean
    row_means          <- rowMeans(expr_mat, na.rm = TRUE)
    expr_mat           <- expr_mat[order(-row_means), , drop = FALSE]
    symbols_ordered    <- symbols[order(-row_means)]
    expr_mat           <- expr_mat[!duplicated(symbols_ordered), , drop = FALSE]
    rownames(expr_mat) <- symbols_ordered[!duplicated(symbols_ordered)]
  }
  
  # Log-transform if data looks like raw intensities (median > 100)
  if (median(expr_mat, na.rm = TRUE) > 100) {
    expr_mat <- log2(expr_mat + 1)
    cat("  Log2-transformed expression (detected raw intensities)\n")
  }
  
  cat(" ", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")
  
  # ── Clinical / survival data ───────────────────────────────────────────────
  pdata <- pData(gse)
  
  # Survival time
  raw_time <- suppressWarnings(as.numeric(pdata[[cohort_cfg$surv_time]]))
  os_time  <- if (cohort_cfg$time_unit == "months") raw_time * 30.44 else raw_time
  
  # Event indicator
  raw_event <- tolower(trimws(pdata[[cohort_cfg$surv_event]]))
  os_event  <- as.integer(raw_event == tolower(cohort_cfg$event_level))
  
  # HPV status (optional)
  hpv_status <- if (!is.null(cohort_cfg$hpv_col) &&
                    cohort_cfg$hpv_col %in% colnames(pdata)) {
    hpv_raw <- tolower(trimws(pdata[[cohort_cfg$hpv_col]]))
    dplyr::case_when(
      grepl("pos|yes|\\+|hpv16|hpv18", hpv_raw) ~ "positive",
      grepl("neg|no|\\-",              hpv_raw) ~ "negative",
      TRUE                                       ~ NA_character_
    )
  } else {
    NA_character_
  }
  
  clin <- data.frame(
    submitter_id = rownames(pdata),
    os_time      = os_time,
    os_event     = os_event,
    hpv_status   = hpv_status,
    stringsAsFactors = FALSE
  ) |> dplyr::filter(!is.na(os_time), os_time > 0)
  
  cat("  Clinical rows retained:", nrow(clin), "\n")
  cat("  HPV breakdown:"); print(table(clin$hpv_status, useNA = "always"))
  
  result <- list(expr_mat = expr_mat, clin = clin)
  saveRDS(result, cache_path)
  return(result)
}


# ── Run survival analysis for one cohort, return Cox summary ─────────────────

#' Score gene sets + run Cox on a single cohort (TCGA or GEO).
#' Wraps run_gsva() and run_cox_survival() so both are called identically
#' for every cohort.
#'
#' @param expr_mat    genes x samples log-normalized matrix
#' @param clin        data.frame(submitter_id, os_time, os_event, hpv_status)
#' @param gene_sets   named list of gene symbol vectors (from extract_go_genesets)
#' @param cohort_name string label used in plot titles and file names
#' @param out_path    directory for outputs
#' @return data.frame of Cox results with cohort column appended
run_cohort_survival <- function(expr_mat, clin, gene_sets,
                                cohort_name, out_path) {
  cohort_dir <- file.path(out_path, cohort_name)
  dir.create(cohort_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\n── Cohort:", cohort_name, "──\n")
  
  # GSVA scoring (with per-cohort cache)
  gsva_cache <- file.path(cohort_dir, "gsva_scores.rds")
  if (file.exists(gsva_cache)) {
    scores_df <- readRDS(gsva_cache)
  } else {
    scores_df <- run_gsva(expr_mat, gene_sets)   # existing function, unchanged
    saveRDS(scores_df, gsva_cache)
  }
  
  # KM curves
  run_km_survival(scores_df, clin, cohort_dir,   # existing function, unchanged
                  stratify_hpv = any(!is.na(clin$hpv_status)))
  
  # Cox — existing function returns invisibly; we capture and tag it
  cox_res <- run_cox_survival(scores_df, clin, cohort_dir)  # existing, unchanged
  cox_res$cohort <- cohort_name
  
  return(cox_res)
}


# ── Forest plot across cohorts ────────────────────────────────────────────────

#' Build a cross-cohort forest plot for each pathway.
#'
#' @param all_cox  data.frame — rbind of run_cohort_survival() outputs,
#'                 must have columns: pathway, HR, CI_lower, CI_upper,
#'                 p_adj, cohort
#' @param out_path directory for PNG / CSV outputs
run_multi_cohort_forest <- function(all_cox, out_path) {
  cat("\nBuilding cross-cohort forest plots...\n")
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(all_cox, file.path(out_path, "all_cohorts_cox.csv"), row.names = FALSE)
  
  pathways <- unique(all_cox$pathway)
  
  for (pw in pathways) {
    df <- all_cox |>
      dplyr::filter(pathway == pw) |>
      dplyr::mutate(
        sig   = dplyr::case_when(
          p_adj < 0.05 & HR > 1 ~ "Risk (FDR<0.05)",
          p_adj < 0.05 & HR < 1 ~ "Protective (FDR<0.05)",
          TRUE                   ~ "NS"
        ),
        # Replicated = same HR direction as TCGA and p_adj < 0.05
        label = paste0(cohort,
                       ifelse(p_adj < 0.05, " *", ""),
                       "  HR=", round(HR, 2))
      )
    
    # Reference HR direction from discovery cohort
    tcga_hr <- df$HR[df$cohort == "TCGA-HNSC"]
    discovery_direction <- if (length(tcga_hr) > 0 && tcga_hr > 1) "risk" else "protective"
    
    p <- ggplot(df, aes(x = HR, y = reorder(cohort, HR),
                        xmin = CI_lower, xmax = CI_upper,
                        color = sig)) +
      geom_pointrange(size = 0.7, fatten = 3) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
      scale_color_manual(values = c(
        "Risk (FDR<0.05)"        = "#E64B35",
        "Protective (FDR<0.05)" = "#4DBBD5",
        "NS"                     = "grey60"
      )) +
      # Annotation: how many cohorts replicate the discovery direction
      annotate("text",
               x     = max(df$CI_upper, na.rm = TRUE) * 0.95,
               y     = 0.6,
               label = paste0("Replication: ",
                              sum(df$cohort != "TCGA-HNSC" &
                                    df$p_adj < 0.05 &
                                    ((discovery_direction == "risk" & df$HR > 1) |
                                       (discovery_direction == "protective" & df$HR < 1))),
                              "/", sum(df$cohort != "TCGA-HNSC"), " cohorts"),
               hjust = 1, size = 3, color = "grey30") +
      labs(
        title  = paste0("Cross-cohort validation — ", gsub("_", " ", pw)),
        x      = "Hazard Ratio (95% CI)",
        y      = NULL,
        color  = NULL,
        caption = "* FDR < 0.05 within cohort"
      ) +
      theme_classic() +
      theme(legend.position = "bottom")
    
    ggsave(
      file.path(out_path, paste0("forest_multicohort_", pw, ".png")),
      p, width = 8, height = max(3, nrow(df) * 0.6 + 2), dpi = 150
    )
  }
  cat("Forest plots saved to:", out_path, "\n")
}

# ── Main entry point ──────────────────────────────────────────────────────────

#' Run the full TCGA survival association pipeline.
#'
#' Call this after hpv_analysis() in st_master.R, passing the go_results object.
#'
#' @param go_results     output of run_go()
#' @param cache_folder   shared cache directory
#' @param picture_folder output directory for plots/CSV
#' @param stratify_hpv   logical — produce HPV+ / HPV- stratified KM plots
#' @param clinical_csv   optional path to pre-processed clinical CSV
#'                       (survival_data_HNSC.csv). When provided, skips the
#'                       TCGAbiolinks clinical download and uses this file
#'                       instead — recommended because it has clean, pre-validated
#'                       HPV status (HNSC_HPV+ / HNSC_HPV-) and pre-computed
#'                       survival times. The expression matrix is still
#'                       downloaded via TCGAbiolinks regardless.
run_tcga_survival <- function(go_results, cache_folder, picture_folder,
                               stratify_hpv = TRUE, clinical_csv = NULL) {

  cat("\n\n::::::::::STARTING TCGA SURVIVAL ASSOCIATION::::::::::\n")
  # ── Download expression data (always via TCGAbiolinks) ────────────────
  tcga     <- get_tcga_hnsc(cache_folder)
  expr_mat <- build_expr_matrix(tcga$expr)

  # ── Clinical / survival data ──────────────────────────────────────────────
  # Prefer the pre-processed CSV when available — it has validated HPV labels
  # and pre-computed survival times, avoiding the ambiguous column parsing
  # in tidy_clinical().
  if (!is.null(clinical_csv) && file.exists(clinical_csv)) {
    cat("Using pre-processed clinical CSV for survival data.\n")
    clin <- load_external_clinical(clinical_csv)
  } else {
    if (!is.null(clinical_csv)) {
      clin <- tidy_clinical(tcga$clinical, subtype = tcga$subtype)
      cat("Warning: clinical_csv path not found —", clinical_csv,
          "\nFalling back to TCGAbiolinks clinical data.\n")
    }
    
  }

  # ── Extract gene sets from spatial GO results ─────────────────────────────
  gene_sets  <- extract_go_genesets(go_results, n_terms = N_GO_TERMS)

  # # ── Score pathways in TCGA ────────────────────────────────────────────────
  # gsva_cache <- file.path(cache_folder, "tcga_gsva_scores.rds")
  # if (file.exists(gsva_cache)) {
  #   cat("Loading GSVA scores from cache...\n")
  #   scores_df <- readRDS(gsva_cache)
  # } else {
  #   scores_df <- run_gsva(expr_mat, gene_sets)
  #   saveRDS(scores_df, gsva_cache)
  #   cat("GSVA scores cached.\n")
  # }
  # 
  # # ── Survival analysis ─────────────────────────────────────────────────────
  # run_km_survival(scores_df, clin, picture_folder, stratify_hpv = stratify_hpv)
  # run_cox_survival(scores_df, clin, picture_folder)

  # ── Discovery cohort ──────────────────────────────────────────────────────
  all_cox <- run_cohort_survival(expr_mat, clin, gene_sets,
                                 cohort_name = "TCGA-HNSC",
                                 out_path    = picture_folder)
  
  # ── External validation ───────────────────────────────────────────────────
  if (validate) {
    cat("\n::::::::::STARTING EXTERNAL VALIDATION::::::::::\n")
    for (cohort_name in names(VALIDATION_COHORTS)) {
      cohort_cfg <- VALIDATION_COHORTS[[cohort_name]]
      geo_data   <- tryCatch(
        load_geo_cohort(cohort_cfg, cache_folder),
        error = function(e) {
          cat("  Skipping", cohort_name, "— download failed:", conditionMessage(e), "\n")
          NULL
        }
      )
      if (is.null(geo_data)) next
      
      cox_res <- run_cohort_survival(
        expr_mat    = geo_data$expr_mat,
        clin        = geo_data$clin,
        gene_sets   = gene_sets,
        cohort_name = cohort_name,
        out_path    = picture_folder
      )
      all_cox <- dplyr::bind_rows(all_cox, cox_res)
    }
    
    # ── Cross-cohort forest plots ─────────────────────────────────────────
    run_multi_cohort_forest(all_cox,
                            out_path = file.path(picture_folder, "multicohort_forest"))
  }  
  
  cat("\n\n::::::::::TCGA SURVIVAL ASSOCIATION COMPLETE::::::::::\n")
}
