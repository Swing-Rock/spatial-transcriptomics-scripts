# ============================================================
# tcga_survival.R
# Associate HPV pathway gene sets (from spatial GO/KEGG enrichment)
# with TCGA-HNSC bulk RNA-seq survival outcomes.
#
# Inputs:
#   go_results  — output of run_go()  from st_hpv_pathway_analysis.R
#   kegg_results— output of run_kegg() (optional, used for gene set labels)
#   cache_folder— directory for .rds checkpoints (shared with main pipeline)
#   picture_folder — directory for PNG/CSV outputs
#
# Workflow:
#   1. Download TCGA-HNSC RNA-seq counts + clinical data  (cached)
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


# ── Main entry point ──────────────────────────────────────────────────────────

#' Run the full TCGA survival association pipeline.
#'
#' Call this after hpv_analysis() in st_master.R, passing the go_results object.
#'
#' @param go_results   output of run_go()
#' @param cache_folder shared cache directory
#' @param picture_folder output directory for plots/CSV
#' @param stratify_hpv   logical — produce HPV+ / HPV- stratified KM plots
run_tcga_survival <- function(go_results, cache_folder, picture_folder,
                               stratify_hpv = TRUE) {

  cat("\n\n::::::::::STARTING TCGA SURVIVAL ASSOCIATION::::::::::\n")

  # ── Download / load data ──────────────────────────────────────────────────
  tcga       <- get_tcga_hnsc(cache_folder)
  expr_mat   <- build_expr_matrix(tcga$expr)
  clin       <- tidy_clinical(tcga$clinical, subtype = tcga$subtype)

  # ── Extract gene sets from spatial GO results ─────────────────────────────
  gene_sets  <- extract_go_genesets(go_results, n_terms = N_GO_TERMS)

  # ── Score pathways in TCGA ────────────────────────────────────────────────
  gsva_cache <- file.path(cache_folder, "tcga_gsva_scores.rds")
  if (file.exists(gsva_cache)) {
    cat("Loading GSVA scores from cache...\n")
    scores_df <- readRDS(gsva_cache)
  } else {
    scores_df <- run_gsva(expr_mat, gene_sets)
    saveRDS(scores_df, gsva_cache)
    cat("GSVA scores cached.\n")
  }

  # ── Survival analysis ─────────────────────────────────────────────────────
  run_km_survival(scores_df, clin, picture_folder, stratify_hpv = stratify_hpv)
  run_cox_survival(scores_df, clin, picture_folder)

  cat("\n\n::::::::::TCGA SURVIVAL ASSOCIATION COMPLETE::::::::::\n")
}
