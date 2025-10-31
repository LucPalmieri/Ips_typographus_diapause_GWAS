#!/usr/bin/env Rscript

################################################################################
# Bayes Factor Manhattan Plot Script
# Purpose: Visualize Bayes Factor values for SNPs across the genome
# Input: Tab-delimited file with CHROM, POS, and BF_median columns
# Output: Manhattan-style plot with facets by chromosome
################################################################################

library(dplyr)
library(ggplot2)

#===============================================================================
# FUNCTION: Create BF Manhattan plot
#===============================================================================
plot_bf_manhattan <- function(input_file,
                              output_file = "BF_manhattan_plot.pdf",
                              bf_column = "BF_median",
                              threshold_moderate = 2,
                              threshold_strong = 6,
                              y_max = 20,
                              plot_width = 12,
                              plot_height = 8,
                              facet_ncol = 6) {
  
  cat("Reading BF data from:", input_file, "\n")
  
  # Read and validate input data
  df <- read.delim(input_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Check required columns
  required_cols <- c("CHROM", "POS", bf_column)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  cat("Processing", nrow(df), "SNPs across", length(unique(df$CHROM)), "chromosomes\n")
  
  # Clean and prepare data
  df$POS <- as.numeric(gsub(",", "", df$POS))
  df[[bf_column]] <- as.numeric(df[[bf_column]])
  
  # Remove rows with missing data
  df <- df[!is.na(df$POS) & !is.na(df[[bf_column]]), ]
  cat("After removing missing data:", nrow(df), "SNPs retained\n")
  
  # Force negative or zero BF values to 0
  df <- df %>%
    mutate(BF_plot = ifelse(.data[[bf_column]] <= 0, 0, .data[[bf_column]]))
  
  # Create faceting variable: keep linkage group (LG) chromosomes, group others as "Uncharted"
  df <- df %>% 
    mutate(facet = ifelse(grepl("^LG", CHROM), CHROM, "Uncharted"))
  
  # Order facets: sort LG chromosomes numerically, then "Uncharted" last
  lg_chrom <- unique(df$CHROM[grepl("^LG", df$CHROM)])
  lg_sorted <- sort(lg_chrom)
  facet_levels <- c(lg_sorted, "Uncharted")
  df$facet <- factor(df$facet, levels = facet_levels)
  
  # Count significant SNPs
  n_moderate <- sum(df$BF_plot >= threshold_moderate & df$BF_plot < threshold_strong)
  n_strong <- sum(df$BF_plot >= threshold_strong)
  cat("\nSummary:\n")
  cat("  SNPs with moderate evidence (BF >=", threshold_moderate, "):", n_moderate, "\n")
  cat("  SNPs with strong evidence (BF >=", threshold_strong, "):", n_strong, "\n")
  
  # Create Manhattan plot
  cat("\nGenerating Manhattan plot...\n")
  p <- ggplot(df, aes(x = POS, y = BF_plot)) +
    geom_point(alpha = 0.7, size = 1.5, color = "black") +
    scale_y_continuous(
      limits = c(0, y_max), 
      breaks = seq(0, y_max, by = 2)
    ) +
    scale_x_continuous(
      breaks = function(x) max(x),
      labels = function(x) paste0(format(x, scientific = FALSE, big.mark = ","), " bp")
    ) +
    geom_hline(
      yintercept = threshold_moderate, 
      color = "red", 
      linetype = "dashed",
      linewidth = 0.5
    ) +
    geom_hline(
      yintercept = threshold_strong, 
      color = "blue", 
      linetype = "dashed",
      linewidth = 0.5
    ) +
    labs(
      x = "Position",
      y = "Bayes Factor (median)",
      title = "Bayes Factor by Position across Chromosomes",
      subtitle = paste0("Red line: BF = ", threshold_moderate, 
                       " (moderate evidence) | Blue line: BF = ", 
                       threshold_strong, " (strong evidence)")
    ) +
    facet_wrap(~facet, scales = "free_x", ncol = facet_ncol) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10)
    )
  
  # Save plot
  ggsave(output_file, plot = p, width = plot_width, height = plot_height)
  cat("Plot saved to:", output_file, "\n")
  
  # Return summary statistics
  return(invisible(list(
    total_snps = nrow(df),
    n_moderate = n_moderate,
    n_strong = n_strong,
    chromosomes = unique(df$CHROM)
  )))
}

#===============================================================================
# FUNCTION: Export significant SNPs to file
#===============================================================================
export_significant_snps <- function(input_file,
                                    bf_column = "BF_median",
                                    threshold = 6,
                                    output_file = "significant_BF_SNPs.txt") {
  
  cat("\nExporting significant SNPs (BF >=", threshold, ")...\n")
  
  # Read data
  df <- read.delim(input_file, header = TRUE, stringsAsFactors = FALSE)
  df$POS <- as.numeric(gsub(",", "", df$POS))
  df[[bf_column]] <- as.numeric(df[[bf_column]])
  
  # Filter significant SNPs
  significant <- df[df[[bf_column]] >= threshold, ]
  significant <- significant[order(-significant[[bf_column]]), ]
  
  cat("Found", nrow(significant), "significant SNPs\n")
  
  # Save to file
  write.table(
    significant,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  cat("Significant SNPs saved to:", output_file, "\n")
  
  return(invisible(significant))
}

#===============================================================================
# MAIN EXECUTION
#===============================================================================

# Example usage
if (sys.nframe() == 0) {
  
  cat("=== BF MANHATTAN PLOT GENERATION ===\n\n")
  
  # Plot with default settings
  results <- plot_bf_manhattan(
    input_file = "BF_position.txt",
    output_file = "BF_manhattan_plot.pdf",
    bf_column = "BF_median",
    threshold_moderate = 2,
    threshold_strong = 6,
    y_max = 20,
    plot_width = 12,
    plot_height = 8,
    facet_ncol = 6
  )
  
  # Export significant SNPs
  significant_snps <- export_significant_snps(
    input_file = "BF_position.txt",
    bf_column = "BF_median",
    threshold = 6,
    output_file = "significant_BF_SNPs.txt"
  )
  
  cat("\nAnalysis complete!\n")
}
