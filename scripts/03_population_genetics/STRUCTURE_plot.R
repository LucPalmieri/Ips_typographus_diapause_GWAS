#!/usr/bin/env Rscript

################################################################################
# STRUCTURE Results Visualization Script
# Purpose: Create publication-quality STRUCTURE bar plots using pophelper
# Input: STRUCTURE output files and population labels
# Output: PDF with ancestry proportion bar plots
################################################################################

library(pophelper)

#===============================================================================
# FUNCTION: Plot STRUCTURE results
#===============================================================================
plot_structure_results <- function(structure_files,
                                   pop_labels_file = NULL,
                                   output_file = "STRUCTURE_plot.pdf",
                                   k_colors = c("#74add1", "#fee090", "#d73027"),
                                   sort_by = "Cluster1",
                                   plot_width = 42,
                                   plot_height = 10,
                                   dpi = 300,
                                   show_dividers = FALSE,
                                   show_legend = TRUE,
                                   show_group_labels = TRUE,
                                   ind_label_angle = 90,
                                   ind_label_size = 3) {
  
  cat("Reading STRUCTURE results from:", structure_files, "\n")
  
  # Read Q-matrix files
  results <- readQ(structure_files)
  
  # Read population labels if provided
  if (!is.null(pop_labels_file)) {
    cat("Reading population labels from:", pop_labels_file, "\n")
    pops <- read.table(pop_labels_file, header = FALSE, stringsAsFactors = FALSE)$V1
    pops_df <- data.frame(Population = pops)
  } else {
    pops_df <- NULL
  }
  
  # Validate K (number of clusters)
  k_value <- ncol(results[[1]])
  cat("Detected K =", k_value, "clusters\n")
  
  # Adjust colors if needed
  if (length(k_colors) < k_value) {
    cat("Warning: Provided", length(k_colors), "colors but K =", k_value, "\n")
    cat("Using default color palette\n")
    k_colors <- NULL
  } else if (length(k_colors) > k_value) {
    cat("Note: Using first", k_value, "colors from provided palette\n")
    k_colors <- k_colors[1:k_value]
  }
  
  # Create STRUCTURE plot
  cat("Generating STRUCTURE plot...\n")
  plotQ(
    results,
    imgtype = "pdf",
    imgoutput = "sep",
    outputfilename = tools::file_path_sans_ext(output_file),
    width = plot_width,
    height = plot_height,
    dpi = dpi,
    barsize = 1,
    showdiv = show_dividers,
    showgrplab = show_group_labels,
    showlegend = show_legend,
    showyaxis = TRUE,
    sortind = sort_by,
    grplab = pops_df,
    ordergrp = TRUE,
    clustercol = k_colors,
    splabsize = 2,
    showindlab = TRUE,
    useindlab = TRUE,
    indlabangle = ind_label_angle,
    indlabsize = ind_label_size,
    indlabheight = 0.3,
    basesize = 11,
    panelspacer = 0.1,
    barbordersize = 0.1,
    barbordercol = "white",
    exportpath = getwd()
  )
  
  cat("Plot saved to:", output_file, "\n")
  
  # Calculate and report ancestry proportions by population
  if (!is.null(pops_df)) {
    cat("\nAncestry proportions by population:\n")
    q_matrix <- as.matrix(results[[1]])
    for (pop in unique(pops)) {
      pop_indices <- which(pops == pop)
      pop_means <- colMeans(q_matrix[pop_indices, , drop = FALSE])
      cat("  ", pop, ": ", paste(sprintf("%.3f", pop_means), collapse = " | "), "\n")
    }
  }
  
  return(invisible(results))
}

#===============================================================================
# FUNCTION: Compare multiple K values
#===============================================================================
compare_k_values <- function(structure_files_list,
                             k_values,
                             pop_labels_file = NULL,
                             output_prefix = "STRUCTURE_K",
                             plot_width = 42,
                             plot_height = 10) {
  
  cat("Comparing STRUCTURE results for K =", paste(k_values, collapse = ", "), "\n")
  
  # Read population labels if provided
  if (!is.null(pop_labels_file)) {
    pops <- read.table(pop_labels_file, header = FALSE, stringsAsFactors = FALSE)$V1
    pops_df <- data.frame(Population = pops)
  } else {
    pops_df <- NULL
  }
  
  # Read all Q-matrices
  all_results <- lapply(structure_files_list, readQ)
  
  # Create combined plot
  cat("Generating combined K comparison plot...\n")
  plotQ(
    qlist = unlist(all_results, recursive = FALSE),
    imgtype = "pdf",
    imgoutput = "join",
    outputfilename = output_prefix,
    width = plot_width,
    height = plot_height * length(k_values),
    dpi = 300,
    barsize = 1,
    showdiv = FALSE,
    showgrplab = TRUE,
    showlegend = TRUE,
    showyaxis = TRUE,
    grplab = pops_df,
    ordergrp = TRUE,
    splabsize = 2,
    showindlab = TRUE,
    useindlab = TRUE,
    indlabangle = 90,
    indlabsize = 3,
    indlabheight = 0.3,
    basesize = 11,
    panelspacer = 0.1,
    barbordersize = 0.1,
    barbordercol = "white",
    exportpath = getwd()
  )
  
  cat("Combined plot saved to:", paste0(output_prefix, ".pdf\n"))
  
  return(invisible(all_results))
}

#===============================================================================
# FUNCTION: Calculate admixture statistics
#===============================================================================
calculate_admixture_stats <- function(structure_files,
                                      pop_labels_file = NULL,
                                      threshold = 0.8) {
  
  cat("\nCalculating admixture statistics...\n")
  cat("Pure ancestry threshold:", threshold, "\n\n")
  
  # Read Q-matrix
  results <- readQ(structure_files)
  q_matrix <- as.matrix(results[[1]])
  k_value <- ncol(q_matrix)
  
  # Identify pure vs admixed individuals
  max_ancestry <- apply(q_matrix, 1, max)
  pure_individuals <- max_ancestry >= threshold
  
  cat("Overall statistics:\n")
  cat("  Total individuals:", nrow(q_matrix), "\n")
  cat("  Pure individuals (>= ", threshold * 100, "%): ", sum(pure_individuals), "\n", sep = "")
  cat("  Admixed individuals: ", sum(!pure_individuals), "\n", sep = "")
  
  # If population labels provided, calculate by population
  if (!is.null(pop_labels_file)) {
    pops <- read.table(pop_labels_file, header = FALSE, stringsAsFactors = FALSE)$V1
    
    cat("\nBy population:\n")
    for (pop in unique(pops)) {
      pop_indices <- which(pops == pop)
      n_pure <- sum(pure_individuals[pop_indices])
      n_total <- length(pop_indices)
      cat("  ", pop, ": ", n_pure, "/", n_total, " pure (", 
          round(n_pure/n_total * 100, 1), "%)\n", sep = "")
    }
  }
  
  # Return statistics
  stats <- data.frame(
    Individual = 1:nrow(q_matrix),
    MaxAncestry = max_ancestry,
    AssignedCluster = apply(q_matrix, 1, which.max),
    Pure = pure_individuals
  )
  
  if (!is.null(pop_labels_file)) {
    stats$Population <- pops
  }
  
  return(invisible(stats))
}

#===============================================================================
# MAIN EXECUTION
#===============================================================================

# Example usage
if (sys.nframe() == 0) {
  
  cat("=== STRUCTURE VISUALIZATION ===\n\n")
  
  # Single K plot
  results <- plot_structure_results(
    structure_files = "struct_phenotyped_f",
    pop_labels_file = "pop_labels.txt",
    output_file = "STRUCTURE_K3_plot.pdf",
    k_colors = c("#74add1", "#fee090", "#d73027"),
    sort_by = "Cluster1",
    plot_width = 42,
    plot_height = 10,
    dpi = 300,
    show_dividers = FALSE,
    show_legend = TRUE,
    show_group_labels = TRUE
  )
  
  # Calculate admixture statistics
  admix_stats <- calculate_admixture_stats(
    structure_files = "struct_phenotyped_f",
    pop_labels_file = "pop_labels.txt",
    threshold = 0.8
  )
  
  # Save admixture statistics
  write.table(
    admix_stats,
    "admixture_statistics.txt",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  cat("\n=== ALTERNATIVE: Compare multiple K values ===\n")
  cat("To compare K=2,3,4, uncomment and modify:\n")
  cat('# compare_results <- compare_k_values(\n')
  cat('#   structure_files_list = list("struct_K2", "struct_K3", "struct_K4"),\n')
  cat('#   k_values = c(2, 3, 4),\n')
  cat('#   pop_labels_file = "pop_labels.txt",\n')
  cat('#   output_prefix = "STRUCTURE_comparison"\n')
  cat('# )\n')
  
  cat("\nAnalysis complete!\n")
}
