################################################################################
## generate pairwise genetic distance for species per timepoint
################################################################################
library(seqinr)
library(dplyr)
library(purrr)

base_dir <- "."
out_dir  <- "../genetic_distance/dist_tables/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
timepoints <- c("T1", "T2", "T3")

get_tp_distances <- function(sp, tp, fasta_all, metadata) {
  
  target_samples <- metadata$Sample_ID[metadata$tp == tp]
   
  matching_idx <- which(names(fasta_all) %in% target_samples)
  if (length(matching_idx) == 0 && "tp" %in% names(metadata)) {
    matching_idx <- which(metadata$tp == tp)
  }
  
  n_seqs <- length(matching_idx)
  
  if (n_seqs < 2) {
    return(tibble(
      species   = sp,
      timepoint = tp,
      distance  = numeric(0)
    ))
  }
  
  sub_fasta <- fasta_all[matching_idx]
  sub_aln   <- as.alignment(
    nb  = n_seqs,
    nam = names(sub_fasta),
    seq = unname(sapply(sub_fasta, paste0, collapse = ""))
  )
  
  d <- dist.alignment(sub_aln, matrix = "identity")
  d_vals <- as.vector(d)
    
  tibble(
    species   = sp,
    timepoint = tp,
    distance  = d_vals
  )
}

all_pairwise_results <- list()

for (i in seq_along(sp46)) {
  sp <- sp46[i]
  cat(sprintf("[%d/%d] Processing %s...\n", i, length(sp46), sp))
  
  aln_file  <- file.path(base_dir, paste0("output_", sp), paste0("s__", sp, ".StrainPhlAn3_concatenated.aln"))
  meta_file <- file.path(base_dir, paste0("output_", sp), "tp.txt")
  
  if (!file.exists(aln_file) || !file.exists(meta_file)) {
    warning(sprintf("Missing files for %s, skipping.", sp))
    next
  }
  
  fasta_all <- tryCatch(read.fasta(aln_file, as.string = FALSE), error = function(e) NULL)
  meta_df   <- tryCatch(read.delim(meta_file, stringsAsFactors = FALSE), error = function(e) NULL)
  
  if (is.null(fasta_all) || is.null(meta_df) || length(fasta_all) == 0) {
    warning(sprintf("Failed reading alignment/metadata for %s", sp))
    next
  }
  
  sp_dists <- map_dfr(timepoints, ~get_tp_distances(sp, .x, fasta_all, meta_df))
  
  if (nrow(sp_dists) > 0) {
    all_pairwise_results[[sp]] <- sp_dists
    
    out_sp_file <- file.path(out_dir, paste0("dist_pairwise_", sp, ".txt"))
    write.table(
      sp_dists %>% mutate(formatted = paste(species, tolower(timepoint), distance, sep = "***")) %>% pull(formatted),
      out_sp_file,
      quote = FALSE, col.names = FALSE, row.names = FALSE
    )
  }
}

combined_pairwise <- bind_rows(all_pairwise_results)
write.table(
  combined_pairwise,
  file.path(base_dir, "dist_pairwise_TP_all.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

summary_stats <- combined_pairwise %>%
  group_by(species, timepoint) %>%
  summarise(
    n_pairs       = n(),
    mean_distance = mean(distance, na.rm = TRUE),
    sd_distance   = sd(distance, na.rm = TRUE),
    .groups       = "drop"
  )

summary_file <- file.path(base_dir, "pairwise_genetic_distance_summary.txt")
write.table(summary_stats, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

print(head(summary_stats))
