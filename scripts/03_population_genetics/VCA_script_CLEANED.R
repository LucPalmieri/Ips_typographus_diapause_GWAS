#!/usr/bin/env Rscript

################################################################################
# Variance Component Analysis (VCA) Script
# Purpose: Calculate genomic relationship matrix and variance components from VCF
# Author: Cleaned and corrected version
################################################################################

library(vcfR)
library(adegenet)
library(rrBLUP)
library(dplyr)

#===============================================================================
# FUNCTION: Perform VCA analysis on VCF data
#===============================================================================
perform_vca <- function(vcf_file, 
                        pheno_file = NULL,
                        pop_file = NULL,
                        output_prefix = "VCA",
                        missing_threshold = 0.90,
                        num_top_vcs = 3) {
  
  cat("Reading VCF file:", vcf_file, "\n")
  vcf_raw <- read.vcfR(vcf_file)
  
  # Calculate missing data per sample
  cat("Calculating missing data per sample...\n")
  gt_all <- extract.gt(vcf_raw, return.alleles = FALSE)
  missing_per_indv <- apply(gt_all, 2, function(x) sum(is.na(x)) / length(x))
  
  # Identify samples exceeding missing data threshold
  samples_to_remove <- names(missing_per_indv[missing_per_indv > missing_threshold])
  
  if (length(samples_to_remove) > 0) {
    cat("Removing", length(samples_to_remove), "samples with >", 
        missing_threshold * 100, "% missing data:\n")
    print(samples_to_remove)
    
    # Correct way to filter vcfR object
    samples_to_keep <- setdiff(colnames(vcf_raw@gt)[-1], samples_to_remove)
    vcf_filtered <- vcf_raw
    vcf_filtered@gt <- vcf_filtered@gt[, c("FORMAT", samples_to_keep)]
  } else {
    cat("No samples exceed missing data threshold\n")
    vcf_filtered <- vcf_raw
  }
  
  # Get final sample list
  samples_vcf_kept <- colnames(vcf_filtered@gt)[-1]
  cat("Number of samples retained:", length(samples_vcf_kept), "\n")
  
  # Convert VCF to genlight
  cat("Converting VCF to genlight object...\n")
  gen <- vcfR2genlight(vcf_filtered)
  
  # Create genotype matrix (impute missing with mean allele frequency)
  cat("Creating genotype matrix and computing relationship matrix...\n")
  gen_mat <- tab(gen, NA.method = "mean")
  
  # Compute genomic relationship matrix (GRM)
  GRM <- A.mat(gen_mat)
  
  # Eigen decomposition
  cat("Performing eigen decomposition...\n")
  eig_GRM <- eigen(GRM)
  vc_scores <- eig_GRM$vectors
  var_explained_vc <- (eig_GRM$values / sum(eig_GRM$values)) * 100
  
  # Create VC score data frame
  vc_df <- data.frame(SampleID = rownames(gen_mat), vc_scores)
  numVC <- ncol(vc_scores)
  colnames(vc_df)[2:(numVC + 1)] <- paste0("VC", 1:numVC)
  
  cat("Variance explained by first 10 VCs (%):\n")
  print(round(var_explained_vc[1:min(10, length(var_explained_vc))], 2))
  
  # Process phenotype data if provided
  if (!is.null(pheno_file)) {
    cat("\nProcessing phenotype data from:", pheno_file, "\n")
    pheno_data <- read.table(pheno_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Subset to samples in VCF
    pheno_data <- pheno_data[pheno_data$SampleID %in% samples_vcf_kept, ]
    cat("Samples with both genotype and phenotype data:", nrow(pheno_data), "\n")
    
    # Merge VC scores with phenotype
    merged_vc <- left_join(vc_df, pheno_data, by = "SampleID")
    
    # Residualize phenotype if phenotype column exists
    if ("Phenotype" %in% colnames(merged_vc)) {
      cat("Residualizing phenotype against top", num_top_vcs, "VCs...\n")
      merged_vc$Phenotype <- as.numeric(merged_vc$Phenotype)
      
      # Check for valid phenotype data
      if (sum(!is.na(merged_vc$Phenotype)) > num_top_vcs) {
        formula_str <- paste("Phenotype ~", paste(paste0("VC", 1:num_top_vcs), collapse = " + "))
        lm_vc <- lm(as.formula(formula_str), data = merged_vc)
        merged_vc$Phenotype_adj <- resid(lm_vc)
        cat("Residualized phenotype created\n")
      } else {
        cat("Warning: Not enough valid phenotype data for residualization\n")
      }
    }
    
    # Save merged data
    output_file <- paste0(output_prefix, "_VC_scores_with_phenotype.txt")
    write.table(merged_vc, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("VC scores with phenotype data saved to:", output_file, "\n")
    
  } else {
    merged_vc <- vc_df
    output_file <- paste0(output_prefix, "_VC_scores.txt")
    write.table(merged_vc, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("VC scores saved to:", output_file, "\n")
  }
  
  # Generate plots if population info provided
  if (!is.null(pop_file)) {
    cat("\nGenerating population structure plots from:", pop_file, "\n")
    pop_info <- read.table(pop_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    names(pop_info) <- c("SampleID", "Population")
    
    # Merge with VC data
    plot_data <- left_join(vc_df, pop_info, by = "SampleID")
    
    # Check if population info is available
    if (sum(!is.na(plot_data$Population)) > 0) {
      # Define colors (customize as needed)
      unique_pops <- sort(unique(plot_data$Population[!is.na(plot_data$Population)]))
      if (length(unique_pops) == 3) {
        popColors <- c("#d73027", "#74add1", "#fee090")
      } else {
        popColors <- rainbow(length(unique_pops))
      }
      plot_data$PopColor <- popColors[match(plot_data$Population, unique_pops)]
      
      # Plot VC1 vs VC2
      pdf(paste0(output_prefix, "_VC1_vs_VC2.pdf"), width = 7, height = 7)
      plot(plot_data$VC1, plot_data$VC2,
           xlab = paste0("VC1 (", round(var_explained_vc[1], 2), "%)"),
           ylab = paste0("VC2 (", round(var_explained_vc[2], 2), "%)"),
           main = "Variance Components - VC1 vs VC2",
           pch = 21,
           bg = plot_data$PopColor,
           col = "black",
           cex = 1.2)
      legend("bottomleft", legend = unique_pops, pch = 21, 
             pt.bg = popColors, col = "black", cex = 0.9)
      dev.off()
      cat("Plot saved to:", paste0(output_prefix, "_VC1_vs_VC2.pdf\n"))
      
      # Multi-panel plot (first 9 VC pairs)
      pdf(paste0(output_prefix, "_multi_panel_VCs.pdf"), width = 12, height = 12)
      par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
      for (i in 1:min(9, numVC - 1)) {
        plot(plot_data[[paste0("VC", i)]],
             plot_data[[paste0("VC", i + 1)]],
             xlab = paste0("VC", i, " (", round(var_explained_vc[i], 2), "%)"),
             ylab = paste0("VC", i + 1, " (", round(var_explained_vc[i + 1], 2), "%)"),
             main = paste0("VC", i, " vs VC", i + 1),
             pch = 21,
             bg = plot_data$PopColor,
             col = "black",
             cex = 0.8)
        if (i == 1) {
          legend("bottomleft", legend = unique_pops, pch = 21, 
                 pt.bg = popColors, col = "black", cex = 0.7)
        }
      }
      dev.off()
      cat("Multi-panel plot saved to:", paste0(output_prefix, "_multi_panel_VCs.pdf\n"))
    } else {
      cat("Warning: No valid population information found in", pop_file, "\n")
    }
  }
  
  cat("\nVCA analysis complete!\n")
  return(invisible(list(
    vc_scores = merged_vc,
    var_explained = var_explained_vc,
    grm = GRM
  )))
}

#===============================================================================
# MAIN EXECUTION
#===============================================================================

# Example usage for phenotyped samples
cat("=== PHENOTYPED SAMPLES VCA ===\n")
pheno_results <- perform_vca(
  vcf_file = "PHENOTYPED.vcf",
  pheno_file = "phenotype.txt",
  pop_file = "phenotyped_3pop.txt",
  output_prefix = "phenotyped",
  missing_threshold = 0.90,
  num_top_vcs = 3
)

# Example usage for wild samples (without phenotype data)
cat("\n=== WILD SAMPLES VCA ===\n")
wild_results <- perform_vca(
  vcf_file = "WILD.recode.vcf",
  pheno_file = NULL,
  pop_file = NULL,
  output_prefix = "wild",
  missing_threshold = 0.90,
  num_top_vcs = 3
)

cat("\nAll analyses completed successfully!\n")
