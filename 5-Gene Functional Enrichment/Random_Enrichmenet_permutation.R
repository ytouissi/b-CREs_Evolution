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

# ── SOURCE ORDER & FILTER ─────────────────────────────────────────────────────
source_order   <- c("GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC", "TF", "WP")

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
cat("  Adult:",  length(real_terms$Adult),  "\n")
cat("  Fetal:",  length(real_terms$Fetal),  "\n")
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
# PART 2 – PERMUTATION FREQUENCY ANALYSIS PER STAGE (batched + resume)
# ══════════════════════════════════════════════════════════════════════════════

n_perm2            <- 1000
batch_size         <- 100
significance_threshold2 <- 0.05
n_cores            <- min(detectCores() - 1, 8)
max_retries        <- 5
initial_backoff    <- 2
output_file        <- "Combined_Enrichment_Results.xlsx"
checkpoint_dir     <- "perm_checkpoints"

if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir)

# ── SOURCE ORDER & FILTER ─────────────────────────────────────────────────────
source_order <- c("GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC", "TF", "WP")

# ── STAGE FILE DEFINITIONS ────────────────────────────────────────────────────
stages <- list(
  Adult = list(
    fg = "adult_HPS-bCREs.txt",
    bg = "adult_non_HPS-bCREs.txt"
  ),
  Fetal = list(
    fg = "fetal_HPS-bCREs.txt",
    bg = "fetal_non_HPS-bCREs.txt"
  ),
  Shared = list(
    fg = "shared_HPS-bCREs.txt",
    bg = "shared_non_HPS_b-CREs.txt"
  )
)

# ── ROTATING QUEUE BUILDER ────────────────────────────────────────────────────
build_rotating_queue <- function(bg_genes, fg_size, n_perms, seed = 42) {
  set.seed(seed)
  draws    <- vector("list", n_perms)
  pool     <- sample(bg_genes)
  pool_pos <- 1L
  for (i in seq_len(n_perms)) {
    draw <- character(0)
    while (length(draw) < fg_size) {
      needed    <- fg_size - length(draw)
      available <- length(pool) - pool_pos + 1L
      if (available >= needed) {
        draw     <- c(draw, pool[pool_pos:(pool_pos + needed - 1L)])
        pool_pos <- pool_pos + needed
      } else {
        draw     <- c(draw, pool[pool_pos:length(pool)])
        pool     <- sample(bg_genes)
        pool_pos <- 1L
      }
    }
    draws[[i]] <- draw
  }
  draws
}

# ── PERMUTATION WORKER ────────────────────────────────────────────────────────
run_perm2 <- function(random_fg, full_bg, organism, threshold, max_retries, initial_backoff) {
  library(gprofiler2)
  perm_res <- NULL
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
  if (!is.null(perm_res) && !is.null(perm_res$result) && nrow(perm_res$result) > 0) {
    r <- perm_res$result
    return(data.frame(
      term_id   = r$term_id,
      term_name = r$term_name,
      source    = r$source,
      term_size = r$term_size,
      p_value   = r$p_value,
      stringsAsFactors = FALSE
    ))
  }
  return(NULL)
}

# ── AGGREGATE HELPER ──────────────────────────────────────────────────────────
aggregate_hits <- function(all_hits) {
  if (is.null(all_hits) || nrow(all_hits) == 0) {
    return(data.frame(
      source            = character(0),
      term_name         = character(0),
      term_id           = character(0),
      term_size         = integer(0),
      times_significant = integer(0),
      best_adj_p_value  = numeric(0),
      stringsAsFactors  = FALSE
    ))
  }
  term_ids <- unique(all_hits$term_id)
  agg <- do.call(rbind, lapply(term_ids, function(tid) {
    rows <- all_hits[all_hits$term_id == tid, ]
    data.frame(
      source            = rows$source[1],
      term_name         = rows$term_name[1],
      term_id           = tid,
      term_size         = rows$term_size[1],
      times_significant = nrow(rows),
      best_adj_p_value  = min(rows$p_value),
      stringsAsFactors  = FALSE
    )
  }))
  
  # Filter to allowed sources only
  agg <- agg[agg$source %in% source_order, ]
  
  # Sort by source order then times_significant desc then best_adj_p_value
  agg$source <- factor(agg$source, levels = source_order)
  agg <- agg[order(agg$source, -agg$times_significant, agg$best_adj_p_value), ]
  agg$source <- as.character(agg$source)
  rownames(agg) <- NULL
  agg
}

# ── VIRIDIS COLOR HELPER ──────────────────────────────────────────────────────
pval_to_hex <- function(pval, x_max) {
  pos <- min(1, max(0, -log10(pval) / x_max))
  col <- viridis(100, direction = -1)[max(1, round(pos * 99) + 1)]
  substr(col, 2, 7)
}

# ── EXCEL WRITER ──────────────────────────────────────────────────────────────
sheet_name_map <- list(
  Adult  = "S8 Random Adult Genes Set",
  Fetal  = "S8 Random Fetal Genes Set",
  Shared = "S8 Random Shared Genes Set"
)

write_excel_combined <- function(perm_df, stage_results, output_file, perms_done) {
  wb <- createWorkbook()
  
  # Global color scale across all p-value columns
  all_pvals <- c(
    perm_df$real_p_adj[is.finite(perm_df$real_p_adj) & perm_df$real_p_adj > 0],
    unlist(lapply(stage_results, function(df) {
      if (is.null(df) || nrow(df) == 0) return(NULL)
      vals <- df$best_adj_p_value
      vals[is.finite(vals) & vals > 0]
    }))
  )
  x_max <- if (length(all_pvals) > 0) -log10(min(all_pvals)) else 1
  
  apply_color <- function(df, sheet_name, col_name) {
    col_idx <- which(colnames(df) == col_name)
    for (row_i in seq_len(nrow(df))) {
      pval <- df[row_i, col_idx]
      if (is.numeric(pval) && is.finite(pval) && pval > 0) {
        hex <- pval_to_hex(pval, x_max)
        sty <- createStyle(fgFill = paste0("#", hex),
                           fontColour = "#FFFFFF", textDecoration = "bold")
        addStyle(wb, sheet_name, sty, rows = row_i + 1, cols = col_idx, gridExpand = FALSE)
      }
    }
  }
  
  # Sheet 1: unchanged permutation output from Part 1
  addWorksheet(wb, "Table S8 Non-HPS Permutation")
  writeData(wb, "Table S8 Non-HPS Permutation", perm_df, colNames = TRUE)
  apply_color(perm_df, "Table S8 Non-HPS Permutation", "real_p_adj")
  
  # Sheets 2-4: permutation frequency results
  for (stage_name in names(stage_results)) {
    sname <- sheet_name_map[[stage_name]]
    df    <- stage_results[[stage_name]]
    if (nrow(df) > 0) df$out_of_permutations <- perms_done
    
    addWorksheet(wb, sname)
    writeData(wb, sname, df, colNames = TRUE)
    if (nrow(df) > 0) apply_color(df, sname, "best_adj_p_value")
  }
  
  saveWorkbook(wb, output_file, overwrite = TRUE)
}

# ── LOAD STAGE INPUTS ONCE ────────────────────────────────────────────────────
stage_inputs <- lapply(stages, function(s) {
  fg_genes <- trimws(readLines(s$fg, warn = FALSE))
  bg_genes <- trimws(readLines(s$bg, warn = FALSE))
  fg_genes <- fg_genes[fg_genes != ""]
  bg_genes <- bg_genes[bg_genes != ""]
  full_bg  <- unique(c(fg_genes, bg_genes))
  list(fg_genes = fg_genes, bg_genes = bg_genes, full_bg = full_bg)
})

for (stage_name in names(stage_inputs)) {
  s <- stage_inputs[[stage_name]]
  cat(stage_name, "-> Foreground:", length(s$fg_genes),
      "| Background:", length(s$full_bg), "\n")
}

# ── PRE-BUILD ROTATING QUEUES ─────────────────────────────────────────────────
cat("\nBuilding rotating gene queues...\n")
rotating_queues <- lapply(names(stages), function(stage_name) {
  s     <- stage_inputs[[stage_name]]
  queue <- build_rotating_queue(s$bg_genes, length(s$fg_genes), n_perm2, seed = 42)
  cat("  [", stage_name, "] Queue ready:", length(queue), "draws of",
      length(s$fg_genes), "genes\n")
  queue
})
names(rotating_queues) <- names(stages)

# ── RESUME ────────────────────────────────────────────────────────────────────
n_batches <- ceiling(n_perm2 / batch_size)

accumulated_hits <- lapply(names(stages), function(stage_name) {
  rds_path <- file.path(checkpoint_dir, paste0(stage_name, "_hits.rds"))
  if (file.exists(rds_path)) {
    cat("  [RESUME] Loading saved hits for", stage_name, "\n")
    readRDS(rds_path)
  } else NULL
})
names(accumulated_hits) <- names(stages)

progress_file <- file.path(checkpoint_dir, "progress.rds")
if (file.exists(progress_file)) {
  perms_done <- readRDS(progress_file)
  if (perms_done >= n_perm2) {
    cat("\nAll", n_perm2, "permutations already complete.\n")
    cat("Final results are in:", output_file, "\n")
    cat("Delete", checkpoint_dir, "to start a fresh run.\n")
    stop("Already complete.", call. = FALSE)
  }
  start_batch <- (perms_done %/% batch_size) + 1
  cat("\nResuming from permutation", perms_done + 1,
      "(batch", start_batch, "/", n_batches, ")\n")
} else {
  perms_done  <- 0
  start_batch <- 1
  cat("\nStarting fresh.\n")
}

# ── BATCH LOOP ────────────────────────────────────────────────────────────────
for (batch in start_batch:n_batches) {
  
  this_batch_size <- min(batch_size, n_perm2 - perms_done)
  batch_start_idx <- perms_done + 1
  perms_done      <- perms_done + this_batch_size
  
  cat("\n══════════════════════════════════════════\n")
  cat("BATCH", batch, "/", n_batches,
      "| Permutations", batch_start_idx, "to", perms_done, "\n")
  cat("══════════════════════════════════════════\n")
  
  for (stage_name in names(stages)) {
    cat("\n  [", stage_name, "] Running", this_batch_size, "permutations...\n")
    
    s           <- stage_inputs[[stage_name]]
    batch_draws <- rotating_queues[[stage_name]][batch_start_idx:perms_done]
    
    cl <- makeCluster(n_cores)
    batch_list <- parLapply(cl, batch_draws, run_perm2,
                            full_bg         = s$full_bg,
                            organism        = organism,
                            threshold       = significance_threshold2,
                            max_retries     = max_retries,
                            initial_backoff = initial_backoff)
    stopCluster(cl)
    
    new_hits <- do.call(rbind, Filter(Negate(is.null), batch_list))
    if (!is.null(new_hits) && nrow(new_hits) > 0) {
      accumulated_hits[[stage_name]] <- rbind(accumulated_hits[[stage_name]], new_hits)
    }
    
    rds_path <- file.path(checkpoint_dir, paste0(stage_name, "_hits.rds"))
    saveRDS(accumulated_hits[[stage_name]], rds_path)
    
    cat("  [", stage_name, "] Done. Total raw hits so far:",
        ifelse(is.null(accumulated_hits[[stage_name]]), 0,
               nrow(accumulated_hits[[stage_name]])), "\n")
  }
  
  cat("\n  Recalculating aggregates and saving Excel...\n")
  
  current_results <- lapply(names(stages), function(stage_name) {
    agg <- aggregate_hits(accumulated_hits[[stage_name]])
    cat("  [", stage_name, "] Unique terms after filtering:", nrow(agg), "\n")
    agg
  })
  names(current_results) <- names(stages)
  
  write_excel_combined(perm_df, current_results, output_file, perms_done)
  saveRDS(perms_done, progress_file)
  cat("  Excel saved ->", output_file,
      "(", perms_done, "/", n_perm2, "permutations complete)\n")
}

cat("\n\nAll", n_perm2, "permutations complete.\n")
cat("Final results saved to:", output_file, "\n")
cat("Checkpoints in:", checkpoint_dir,
    "— delete this folder to start a fresh run.\n")