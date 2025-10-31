#!/usr/bin/env Rscript

################################################################################
# Binomial Test Script for Classification Accuracy
# Purpose: Statistical testing of observed vs expected classification rates
# Use case: Validate classifier performance against expected proportions
################################################################################

#===============================================================================
# FUNCTION: Perform binomial test with detailed output
#===============================================================================
perform_binomial_test <- function(observed_success,
                                  total_trials,
                                  expected_proportion,
                                  test_name = "Test",
                                  alternative = "two.sided",
                                  confidence_level = 0.95) {
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat(test_name, "\n")
  cat(rep("=", 70), "\n", sep = "")
  
  # Perform binomial test
  test_result <- binom.test(
    x = observed_success,
    n = total_trials,
    p = expected_proportion,
    alternative = alternative,
    conf.level = confidence_level
  )
  
  # Calculate observed proportion
  observed_proportion <- observed_success / total_trials
  
  # Print formatted results
  cat("Observed successes:     ", observed_success, "/", total_trials, 
      " (", round(observed_proportion * 100, 2), "%)\n", sep = "")
  cat("Expected proportion:    ", round(expected_proportion * 100, 2), "%\n", sep = "")
  cat("P-value:                ", format.pval(test_result$p.value, digits = 3), "\n", sep = "")
  cat("Confidence interval:    [", 
      round(test_result$conf.int[1], 4), ", ",
      round(test_result$conf.int[2], 4), "]\n", sep = "")
  
  # Interpretation
  alpha <- 1 - confidence_level
  if (test_result$p.value < alpha) {
    cat("Result:                 SIGNIFICANT (p < ", alpha, ")\n", sep = "")
    cat("Interpretation:         Observed proportion differs significantly from expected\n")
  } else {
    cat("Result:                 NOT SIGNIFICANT (p >= ", alpha, ")\n", sep = "")
    cat("Interpretation:         Observed proportion is consistent with expected\n")
  }
  
  cat("\n")
  
  # Return result invisibly
  return(invisible(test_result))
}

#===============================================================================
# FUNCTION: Batch binomial tests for multiple scenarios
#===============================================================================
batch_binomial_tests <- function(test_scenarios,
                                 output_file = "binomial_test_results.txt") {
  
  cat("\n=== BATCH BINOMIAL TESTING ===\n")
  cat("Running", nrow(test_scenarios), "tests...\n")
  
  # Initialize results storage
  results_list <- list()
  
  # Open connection to output file
  sink(output_file)
  cat("Binomial Test Results\n")
  cat("Generated:", date(), "\n\n")
  
  # Run each test
  for (i in 1:nrow(test_scenarios)) {
    scenario <- test_scenarios[i, ]
    
    result <- perform_binomial_test(
      observed_success = scenario$observed_success,
      total_trials = scenario$total_trials,
      expected_proportion = scenario$expected_proportion,
      test_name = scenario$test_name,
      alternative = "two.sided"
    )
    
    results_list[[i]] <- result
  }
  
  # Close output file
  sink()
  
  cat("Results saved to:", output_file, "\n")
  
  return(invisible(results_list))
}

#===============================================================================
# FUNCTION: Create summary table of test results
#===============================================================================
create_summary_table <- function(test_scenarios, test_results) {
  
  summary_df <- data.frame(
    Test = test_scenarios$test_name,
    Observed = test_scenarios$observed_success,
    Total = test_scenarios$total_trials,
    Expected_Prop = test_scenarios$expected_proportion,
    Observed_Prop = test_scenarios$observed_success / test_scenarios$total_trials,
    P_value = sapply(test_results, function(x) x$p.value),
    CI_lower = sapply(test_results, function(x) x$conf.int[1]),
    CI_upper = sapply(test_results, function(x) x$conf.int[2]),
    Significant = sapply(test_results, function(x) x$p.value < 0.05)
  )
  
  return(summary_df)
}

#===============================================================================
# FUNCTION: Visualize test results
#===============================================================================
plot_binomial_results <- function(summary_table,
                                  output_file = "binomial_test_plot.pdf") {
  
  library(ggplot2)
  
  # Prepare data for plotting
  summary_table$Test <- factor(summary_table$Test, levels = summary_table$Test)
  
  # Create plot
  p <- ggplot(summary_table, aes(x = Test, y = Observed_Prop)) +
    geom_point(size = 4, aes(color = Significant)) +
    geom_errorbar(
      aes(ymin = CI_lower, ymax = CI_upper, color = Significant),
      width = 0.2,
      linewidth = 1
    ) +
    geom_hline(
      aes(yintercept = Expected_Prop),
      linetype = "dashed",
      color = "gray50"
    ) +
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "black"),
      labels = c("Not Significant", "Significant (p < 0.05)")
    ) +
    labs(
      x = "Test Scenario",
      y = "Proportion",
      title = "Binomial Test Results: Observed vs Expected Proportions",
      subtitle = "Error bars show 95% confidence intervals; dashed line shows expected proportion",
      color = "Statistical\nSignificance"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )
  
  # Save plot
  ggsave(output_file, plot = p, width = 10, height = 6)
  cat("Plot saved to:", output_file, "\n")
  
  return(invisible(p))
}

#===============================================================================
# MAIN EXECUTION: Mock classifier validation tests
#===============================================================================

if (sys.nframe() == 0) {
  
  cat("\n", rep("#", 70), "\n", sep = "")
  cat("BINOMIAL TESTS: CLASSIFIER VALIDATION\n")
  cat(rep("#", 70), "\n", sep = "")
  cat("\nTesting classifier performance against expected phenotype proportions\n")
  cat("Scenario: Mock datasets with known phenotype ratios (100%, 75%, 50%, 25%, 0%)\n")
  
  # Define test scenarios (from original script)
  test_scenarios <- data.frame(
    test_name = c("Mock100 (Expected: 100% Phenotype)",
                  "Mock75 (Expected: 75% Phenotype)",
                  "Mock50 (Expected: 50% Phenotype)",
                  "Mock25 (Expected: 25% Phenotype)",
                  "Mock0 (Expected: 0% Phenotype)"),
    observed_success = c(97, 74, 49, 27, 2),
    total_trials = c(100, 100, 100, 100, 100),
    expected_proportion = c(1.0, 0.75, 0.50, 0.25, 0.0),
    stringsAsFactors = FALSE
  )
  
  # Run all tests and save to file
  test_results <- batch_binomial_tests(
    test_scenarios = test_scenarios,
    output_file = "binomial_test_results.txt"
  )
  
  # Create summary table
  summary_table <- create_summary_table(test_scenarios, test_results)
  
  # Display summary
  cat("\n=== SUMMARY TABLE ===\n\n")
  print(summary_table, row.names = FALSE, digits = 4)
  
  # Save summary table
  write.table(
    summary_table,
    "binomial_test_summary.txt",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  cat("\nSummary table saved to: binomial_test_summary.txt\n")
  
  # Create visualization
  cat("\nGenerating visualization...\n")
  plot_binomial_results(
    summary_table = summary_table,
    output_file = "binomial_test_plot.pdf"
  )
  
  # Overall assessment
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("OVERALL ASSESSMENT\n")
  cat(rep("=", 70), "\n", sep = "")
  
  n_significant <- sum(summary_table$Significant)
  n_total <- nrow(summary_table)
  
  cat("Tests conducted:        ", n_total, "\n", sep = "")
  cat("Significant deviations: ", n_significant, "\n", sep = "")
  cat("Non-significant:        ", n_total - n_significant, "\n", sep = "")
  
  if (n_significant == 0) {
    cat("\nConclusion: Classifier performs as expected across all test scenarios.\n")
  } else if (n_significant == n_total) {
    cat("\nConclusion: Classifier shows significant bias in all test scenarios.\n")
    cat("            Further investigation recommended.\n")
  } else {
    cat("\nConclusion: Classifier shows mixed performance.\n")
    cat("            Significant deviations in ", n_significant, " of ", n_total, " scenarios.\n", sep = "")
  }
  
  cat("\nAnalysis complete!\n")
}
