# ── TOGGLE ───────────────────────────────────────────────────────────────────
USE_EXISTING_RANDOM <- TRUE   # TRUE  = use existing random_*_genes.txt
                               # FALSE = draw new random samples
# ─────────────────────────────────────────────────────────────────────────────

required_packages <- c("gprofiler2", "readxl", "openxlsx", "ggplot2", "parallel", "httr", "viridis")
installed <- rownames(installed.packages())
to_install <- required_packages[!required_packages %in% installed]
if (length(to_install) > 0) {
  message("Installing missing packages: ", paste(to_install, collapse = ", "))
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

library(gprofiler2)
library(readxl)
library(openxlsx)
library(ggplot2)
library(parallel)
library(httr)
library(viridis)

httr::set_config(httr::timeout(120))
set.seed(42)

# ══════════════════════════════════════════════════════════════════════════════
# PART 1 – PERMUTATION SANITY CHECK (from Table_S7)
# ══════════════════════════════════════════════════════════════════════════════

n_permutations         <- 1000
significance_threshold <- 0.05
organism               <- "hsapiens"
n_cores                <- min(detectCores() - 1, 30)
max_retries            <- 5
initial_backoff        <- 2

cat("Using", n_cores, "cores\n")

run_one_permutation <- function(i, bg_genes, fg_size, full_bg, target_terms,
                                organism, threshold, max_retries, initial_backoff) {
  library(gprofiler2)
  random_fg <- sample(bg_genes, size = fg_size, replace = FALSE)
  perm_res  <- NULL
  for (attempt in seq_len(max_retries)) {
    perm_res <- tryCatch({
      gost(query = random_fg, organism = organism, significant = TRUE,
           user_threshold = threshold, correction_method = "g_SCS",
           domain_scope = "custom", custom_bg = full_bg)
    }, error = function(e) NULL)
    if (!is.null(perm_res)) break
    Sys.sleep(initial_backoff * (2 ^ (attempt - 1)) + runif(1, 0, 2))
  }
  Sys.sleep(runif(1, 0.5, 1.5))
  if (!is.null(perm_res) && !is.null(perm_res$result))
    return(target_terms[target_terms %in% perm_res$result$term_id])
  character(0)
}

s7 <- read_excel("Table_S7.xlsx", sheet = "Table S7")

real_terms <- list(
  Adult  = s7$term_id[!is.na(s7$`adjusted_p_value__HPS Adult b-CREs`)   & s7$`adjusted_p_value__HPS Adult b-CREs`   < 0.05],
  Fetal  = s7$term_id[!is.na(s7$`adjusted_p_value__HPS Fetal b-CREs`)   & s7$`adjusted_p_value__HPS Fetal b-CREs`   < 0.05],
  Shared = s7$term_id[!is.na(s7$`adjusted_p_value__HPS Shared  b-CREs`) & s7$`adjusted_p_value__HPS Shared  b-CREs` < 0.05]
)

cat("Significant terms per stage:\n")
cat("  Adult:", length(real_terms$Adult), "\n")
cat("  Fetal:", length(real_terms$Fetal), "\n")
cat("  Shared:", length(real_terms$Shared), "\n\n")

perm_stages <- list(
  Adult  = list(fg = "adult_HPS-bCREs.txt",  bg = "adult_non_HPS-bCREs.txt"),
  Fetal  = list(fg = "fetal_HPS-bCREs.txt",  bg = "fetal_non_HPS-bCREs.txt"),
  Shared = list(fg = "shared_HPS-bCREs.txt", bg = "shared_non_HPS_b-CREs.txt")
)

all_term_results <- list()

for (stage_name in names(perm_stages)) {

  cat("\n===", stage_name, "===\n")

  fg_genes     <- trimws(readLines(perm_stages[[stage_name]]$fg, warn = FALSE))
  bg_genes     <- trimws(readLines(perm_stages[[stage_name]]$bg, warn = FALSE))
  fg_genes     <- fg_genes[fg_genes != ""]
  bg_genes     <- bg_genes[bg_genes != ""]
  full_bg      <- unique(c(fg_genes, bg_genes))
  target_terms <- real_terms[[stage_name]]

  cat("Foreground:", length(fg_genes), "| Background:", length(full_bg),
      "| Tracking:", length(target_terms), "terms\n")
  cat("Running", n_permutations, "permutations across", n_cores, "cores...\n")

  cl <- makeCluster(n_cores)
  clusterSetRNGStream(cl, iseed = 42)
  perm_hits <- parLapply(cl, seq_len(n_permutations), run_one_permutation,
                         bg_genes = bg_genes, fg_size = length(fg_genes),
                         full_bg = full_bg, target_terms = target_terms,
                         organism = organism, threshold = significance_threshold,
                         max_retries = max_retries, initial_backoff = initial_backoff)
  stopCluster(cl)
  cat("  Done.\n")

  term_hit_count <- setNames(rep(0, length(target_terms)), target_terms)
  for (hits in perm_hits) term_hit_count[hits] <- term_hit_count[hits] + 1

  p_col <- switch(stage_name,
    Adult  = "adjusted_p_value__HPS Adult b-CREs",
    Fetal  = "adjusted_p_value__HPS Fetal b-CREs",
    Shared = "adjusted_p_value__HPS Shared  b-CREs"
  )

  stage_results <- data.frame(
    stage           = stage_name,
    source          = s7$source[match(target_terms, s7$term_id)],
    term_id         = target_terms,
    term_name       = s7$term_name[match(target_terms, s7$term_id)],
    real_p_adj      = s7[[p_col]][match(target_terms, s7$term_id)],
    times_in_random = as.numeric(term_hit_count[target_terms]),
    empirical_p     = (as.numeric(term_hit_count[target_terms]) + 1) / (n_permutations + 1),
    stringsAsFactors = FALSE
  )
  stage_results <- stage_results[order(stage_results$empirical_p), ]
  all_term_results[[stage_name]] <- stage_results

  cat("\n  Robust terms (never/rarely in random):\n")
  rare <- stage_results[stage_results$times_in_random <= 5, ]
  if (nrow(rare) > 0) print(head(rare[, c("term_id","term_name","real_p_adj","times_in_random","empirical_p")], 15))

  cat("\n  Non-specific terms (frequent in random):\n")
  common <- stage_results[stage_results$times_in_random >= 50, ]
  if (nrow(common) > 0) print(head(common[, c("term_id","term_name","real_p_adj","times_in_random","empirical_p")], 15))
  else cat("  None\n")
}

perm_df <- do.call(rbind, all_term_results)
write.csv(perm_df, "term_specific_permutation_results.csv", row.names = FALSE)

p <- ggplot(perm_df, aes(x = empirical_p)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ stage, ncol = 1, scales = "free_y") +
  labs(title = "Term-Specific Permutation Test",
       subtitle = "Empirical p-value: fraction of random draws where each term was significant",
       x = "Empirical p-value", y = "Number of terms") +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(face = "bold"))
ggsave("term_specific_permutation_test.pdf", p, width = 8, height = 10)
ggsave("term_specific_permutation_test.png", p, width = 8, height = 10, dpi = 300)

cat("\n========== PERMUTATION SUMMARY ==========\n")
for (stage_name in names(all_term_results)) {
  res <- all_term_results[[stage_name]]
  n_robust <- sum(res$empirical_p < 0.05)
  cat(stage_name, ": ", n_robust, "/", nrow(res),
      " (", round(100 * n_robust / nrow(res), 1), "%) robust\n", sep = "")
}

# ══════════════════════════════════════════════════════════════════════════════
# PART 2 – RANDOM ENRICHMENT PER STAGE
# ══════════════════════════════════════════════════════════════════════════════

rand_stages <- list(
  Adult = list(
    random_fg    = "random_adult_genes.txt",
    non_selected = "adult_non_HPS-bCREs.txt",
    selected     = "adult_HPS-bCREs.txt"
  ),
  Fetal = list(
    random_fg    = "random_fetal_genes.txt",
    non_selected = "fetal_non_HPS-bCREs.txt",
    selected     = "fetal_HPS-bCREs.txt"
  ),
  Shared = list(
    random_fg    = "random_shared_genes.txt",
    non_selected = "shared_non_HPS_b-CREs.txt",
    selected     = "shared_HPS-bCREs.txt"
  )
)

random_results <- list()

for (stage_name in names(rand_stages)) {

  cat("\n===", stage_name, "(random enrichment) ===\n")

  non_selected <- trimws(readLines(rand_stages[[stage_name]]$non_selected, warn = FALSE))
  selected     <- trimws(readLines(rand_stages[[stage_name]]$selected,     warn = FALSE))
  non_selected <- non_selected[non_selected != ""]
  selected     <- selected[selected != ""]

  if (USE_EXISTING_RANDOM) {
    random_fg <- trimws(readLines(rand_stages[[stage_name]]$random_fg, warn = FALSE))
    random_fg <- random_fg[random_fg != ""]
    cat("  Using existing random genes:", length(random_fg), "\n")
  } else {
    pool      <- setdiff(non_selected, selected)
    n         <- min(length(selected), length(pool))
    random_fg <- sample(pool, size = n, replace = FALSE)
    writeLines(random_fg, rand_stages[[stage_name]]$random_fg)
    cat("  New random sample:", length(random_fg), "-> saved to", rand_stages[[stage_name]]$random_fg, "\n")
  }

  full_bg <- unique(c(selected, non_selected))
  cat("  Full bg:", length(full_bg), "\n")

  res <- tryCatch({
    gost(query = random_fg, organism = "hsapiens", significant = TRUE,
         user_threshold = 0.05, correction_method = "g_SCS",
         domain_scope = "custom", custom_bg = full_bg)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(res)) {
    meta <- res$meta$genes_metadata$query$query_1
    cat("  Genes recognised:", meta$ensembl$recognise_numbers,
        "| Not recognised:", meta$ensembl$unrecognise_numbers, "\n")
  }

  label <- paste0("Random ", stage_name, " Non-HPS b-CREs")

  if (!is.null(res) && !is.null(res$result) && nrow(res$result) > 0) {
    df_raw <- res$result[order(res$result$p_value), ]
    df <- data.frame(
      term_name             = df_raw$term_name,
      term_id               = df_raw$term_id,
      term_size             = df_raw$term_size,
      effective_domain_size = df_raw$effective_domain_size,
      adjusted_p_value      = df_raw$p_value,
      query_size            = df_raw$query_size,
      intersection_size     = df_raw$intersection_size,
      stringsAsFactors      = FALSE
    )
    colnames(df)[5] <- paste0("adjusted_p_value__", label)
    colnames(df)[6] <- paste0("query_size__", label)
    colnames(df)[7] <- paste0("intersection_size__", label)
    cat("  Significant terms:", nrow(df), "\n")
    random_results[[stage_name]] <- df
  } else {
    cat("  No significant terms found.\n")
    random_results[[stage_name]] <- NULL
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# PART 3 – WRITE EXCEL WITH VIRIDIS COLORING
# ══════════════════════════════════════════════════════════════════════════════

pval_to_hex <- function(pval, x_max) {
  pos <- min(1, max(0, -log10(pval) / x_max))
  col <- viridis(100, direction = -1)[max(1, round(pos * 99) + 1)]
  substr(col, 2, 7)
}

all_pvals_for_scale <- c(
  perm_df$empirical_p[is.finite(perm_df$empirical_p) & perm_df$empirical_p > 0],
  perm_df$real_p_adj[is.finite(perm_df$real_p_adj)   & perm_df$real_p_adj   > 0],
  unlist(lapply(random_results, function(df) {
    if (is.null(df)) return(NULL)
    pcol <- grep("adjusted_p_value", colnames(df), value = TRUE)
    vals <- df[[pcol]]
    vals[is.finite(vals) & vals > 0]
  }))
)
x_max <- if (length(all_pvals_for_scale) > 0) -log10(min(all_pvals_for_scale)) else 1
cat("\nGlobal color scale: most significant p =", min(all_pvals_for_scale), "\n")

apply_color <- function(wb, sheet_name, df, col_indices) {
  for (ci in col_indices) {
    for (row_i in seq_len(nrow(df))) {
      pval <- df[row_i, ci]
      if (is.numeric(pval) && is.finite(pval) && pval > 0) {
        hex <- pval_to_hex(pval, x_max)
        sty <- createStyle(fgFill = paste0("#", hex),
                           fontColour = "#FFFFFF",
                           textDecoration = "bold")
        addStyle(wb, sheet_name, sty, rows = row_i + 1, cols = ci, gridExpand = FALSE)
      }
    }
  }
}

wb <- createWorkbook()

# Sheet 1: permutation
addWorksheet(wb, "Permutation")
writeData(wb, "Permutation", perm_df)
apply_color(wb, "Permutation", perm_df,
            which(colnames(perm_df) %in% c("real_p_adj", "empirical_p")))

# Sheets 2-4: random enrichment
for (stage_name in names(rand_stages)) {
  df <- random_results[[stage_name]]
  addWorksheet(wb, stage_name)
  if (is.null(df)) {
    writeData(wb, stage_name, data.frame(message = "No significant enrichment terms found"))
  } else {
    writeData(wb, stage_name, df)
    apply_color(wb, stage_name, df, grep("adjusted_p_value", colnames(df)))
  }
}

saveWorkbook(wb, "Combined_Enrichment_Results.xlsx", overwrite = TRUE)
cat("\nDone. Results saved to Combined_Enrichment_Results.xlsx\n")
