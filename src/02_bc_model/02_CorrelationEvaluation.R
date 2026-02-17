
library(grid)
library(ggplot2)
library(reshape2)
library(gridExtra)

dimtest <- function(base_input_dir, 
                    k_list = 2:8, 
                    experiment_range = 2:6,
                    n_clusters = 3,
                    output_subdir = "correlation_results") {
  
  # Validate inputs
  if (!dir.exists(base_input_dir)) {
    stop("Base input directory does not exist: ", base_input_dir)
  }
  
  if (length(k_list) == 0 || length(experiment_range) == 0 || n_clusters <= 0) {
    stop("Invalid parameters: k_list, experiment_range must be non-empty, n_clusters must be positive")
  }
  
  start <- Sys.time()
  cat("Starting correlation analysis...\n")
  cat("Base directory:", base_input_dir, "\n")
  cat("K values:", paste(k_list, collapse = ", "), "\n")
  cat("Experiments:", paste(experiment_range, collapse = ", "), "\n")
  cat("Number of clusters:", n_clusters, "\n\n")
  

  correlation_array <- array(NA, 
                             dim = c(length(experiment_range), length(k_list), n_clusters),
                             dimnames = list(
                               paste0("Exp_", experiment_range),
                               paste0("K_", k_list),
                               paste0("Cluster_", 0:(n_clusters-1))
                             ))
  

  total_files <- length(experiment_range) * length(k_list) * n_clusters
  processed_files <- 0
  successful_correlations <- 0
  
  # Loop through experiments
  for(exp_idx in seq_along(experiment_range)) {
    exp_num <- experiment_range[exp_idx]
    input_dir <- file.path(base_input_dir, as.character(exp_num))
    
    if (!dir.exists(input_dir)) {
      warning("Experiment directory does not exist: ", input_dir)
      next
    }
    
    cat("Processing Experiment", exp_num, "...\n")
    
    for (k_idx in seq_along(k_list)) {
      k <- k_list[k_idx]
      location <- file.path(input_dir, paste0("K=", k), paste0("K=", k))
      
      if (!dir.exists(location)) {
        warning("K directory does not exist: ", location)
        next
      }
      
      for (cluster_id in 0:(n_clusters-1)) {
        processed_files <- processed_files + 1
        
        model_file <- file.path(location, paste0("bc_model_cluster_", cluster_id, ".rds"))
        test_data_file <- file.path(location, paste0("test_data_cluster_", cluster_id, ".rds"))
        
        cat(sprintf("  [%d/%d] Processing K=%d, Cluster %d\n", 
                    processed_files, total_files, k, cluster_id))
        
        # Load model file
        model_data <- tryCatch({
          if (!file.exists(model_file)) {
            stop("Model file does not exist")
          }
          readRDS(model_file)
        }, error = function(e) {
          message("    Error reading model file: ", model_file, " - ", e$message)
          return(NULL)
        })
        
        # Load test data file
        test_data <- tryCatch({
          if (!file.exists(test_data_file)) {
            stop("Test data file does not exist")
          }
          readRDS(test_data_file)
        }, error = function(e) {
          message("    Error reading test file: ", test_data_file, " - ", e$message)
          return(NULL)
        })
        
        # Process correlation if both files loaded successfully
        if (!is.null(model_data) && !is.null(test_data)) {
          tryCatch({
            # Get model matrix
            if (!"UVPM" %in% names(model_data)) {
              stop("UVPM not found in model data")
            }
            model <- model_data$UVPM
            
            # Compute average test matrix
            if (length(test_data) == 0) {
              stop("Test data is empty")
            }
            x_test <- Reduce("+", test_data) / length(test_data)
            
            # Validate matrices
            if (is.null(x_test) || is.null(model)) {
              stop("One or both matrices are NULL")
            }
            
            if (!all(dim(x_test) == dim(model))) {
              stop("Matrix dimensions do not match")
            }
            
            if (!is.matrix(x_test) || !is.matrix(model)) {
              stop("Data is not in matrix format")
            }
            
            # Convert matrices to vectors - only upper triangle to avoid redundancy
            vec1 <- x_test[upper.tri(x_test)]
            vec2 <- model[upper.tri(model)]
            
            if (length(vec1) == 0 || length(vec2) == 0) {
              stop("No upper triangle elements found")
            }
            
          
            correlation_value <- cor(vec1, vec2, use = "complete.obs")
            
            if (is.na(correlation_value)) {
              warning("    Correlation resulted in NA")
            } else {
              correlation_array[exp_idx, k_idx, cluster_id + 1] <- correlation_value
              successful_correlations <- successful_correlations + 1
              cat(sprintf("    Correlation: %.4f\n", correlation_value))
            }
            
          }, error = function(e) {
            message("    Error computing correlation: ", e$message)
          })
        }
        
      
        if (processed_files %% 10 == 0) {
          gc()
        }
      }
    }
  }
  
  cat("\nProcessing Summary:\n")
  cat("Total files processed:", processed_files, "\n")
  cat("Successful correlations:", successful_correlations, "\n")
  cat("Success rate:", sprintf("%.1f%%", 100 * successful_correlations / processed_files), "\n\n")
  
  
  if (successful_correlations == 0) {
    stop("No successful correlations computed. Please check your data and file paths.")
  }
  
  # Calculate means and standard deviations across experiments
  mean_correlations <- apply(correlation_array, c(2, 3), mean, na.rm = TRUE)
  sd_correlations <- apply(correlation_array, c(2, 3), sd, na.rm = TRUE)
  
  # Create summary dataframes
  mean_df <- as.data.frame(mean_correlations)
  colnames(mean_df) <- paste0("Cluster_", 0:(n_clusters-1))
  mean_df$K_value <- k_list
  
  sd_df <- as.data.frame(sd_correlations)
  colnames(sd_df) <- paste0("Cluster_", 0:(n_clusters-1))
  sd_df$K_value <- k_list
  
  # Create output directory
  output_dir <- file.path(base_input_dir, output_subdir)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("Generating visualizations and saving results...\n")
  

  pdf_file <- file.path(output_dir, "Correlation_Analysis_Results.pdf")
  pdf(pdf_file, width = 12, height = 15)
  
  tryCatch({
  
    mean_df_long <- melt(mean_df, id.vars = "K_value", 
                         variable.name = "Cluster", value.name = "Correlation")
    
    p1 <- ggplot(mean_df_long, aes(x = K_value, y = Correlation, color = Cluster)) +
      geom_line(size = 1) +
      geom_point(size = 2) +
      theme_minimal() +
      labs(title = paste("Mean Correlation Coefficients by K-value and Cluster"),
           subtitle = paste("Experiments", min(experiment_range), "-", max(experiment_range)),
           x = "K Value", y = "Mean Correlation") +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 12)) +
      scale_x_continuous(breaks = k_list)
    

    mean_table <- tableGrob(round(mean_df, 4),
                            rows = NULL,
                            cols = c("K", paste("Cluster", 0:(n_clusters-1))),
                            theme = ttheme_minimal(base_size = 10))
    
    sd_table <- tableGrob(round(sd_df, 4),
                          rows = NULL,
                          cols = c("K", paste("Cluster", 0:(n_clusters-1))),
                          theme = ttheme_minimal(base_size = 10))
    
  
    grid.arrange(
      textGrob(paste("Correlation Analysis - Experiments", 
                     min(experiment_range), "to", max(experiment_range)), 
               gp = gpar(fontsize = 16, fontface = "bold")),
      p1,
      textGrob("Mean Correlation Values", 
               gp = gpar(fontsize = 14, fontface = "bold")),
      mean_table,
      textGrob("Standard Deviations", 
               gp = gpar(fontsize = 14, fontface = "bold")),
      sd_table,
      ncol = 1,
      heights = c(0.5, 3, 0.5, 2, 0.5, 2)
    )
    
    
    for(k_idx in seq_along(k_list)) {
      grid.newpage()
      k <- k_list[k_idx]
      k_data <- correlation_array[, k_idx, ]
      colnames(k_data) <- paste0("Cluster_", 0:(n_clusters-1))
      rownames(k_data) <- paste0("Exp_", experiment_range)
      
      grid.text(paste0("Detailed Results for K=", k), 
                y = 0.9, gp = gpar(fontsize = 14, fontface = "bold"))
      
      
      if (any(!is.na(k_data))) {
        grid.table(round(k_data, 4))
      } else {
        grid.text("No valid correlations for this K value", 
                  y = 0.5, gp = gpar(fontsize = 12, col = "red"))
      }
    }
    
  }, error = function(e) {
    message("Error generating plots: ", e$message)
  })
  
  dev.off()
  cat("PDF saved to:", pdf_file, "\n")
  
  # Save raw data
  data_file <- file.path(output_dir, "correlation_analysis_data.rds")
  saveRDS(list(
    correlation_array = correlation_array,
    mean_correlations = mean_correlations,
    sd_correlations = sd_correlations,
    parameters = list(
      k_list = k_list,
      experiment_range = experiment_range,
      n_clusters = n_clusters,
      base_input_dir = base_input_dir
    ),
    summary = list(
      total_files = total_files,
      processed_files = processed_files,
      successful_correlations = successful_correlations,
      success_rate = successful_correlations / processed_files
    )
  ), file = data_file)
  cat("Raw data saved to:", data_file, "\n")
  
 
  write.csv(mean_df, file.path(output_dir, "mean_correlations.csv"), row.names = FALSE)
  write.csv(sd_df, file.path(output_dir, "sd_correlations.csv"), row.names = FALSE)
  

  end <- Sys.time()
  duration <- difftime(end, start, units = "mins")
  cat(sprintf("\nExecution completed in %.2f minutes\n", as.numeric(duration)))
  
  # Print summary statistics
  cat("\nSummary Statistics:\n")
  overall_mean <- mean(correlation_array, na.rm = TRUE)
  overall_sd <- sd(correlation_array, na.rm = TRUE)
  cat("Overall mean correlation:", sprintf("%.4f", overall_mean), "\n")
  cat("Overall standard deviation:", sprintf("%.4f", overall_sd), "\n")
  cat("Range:", sprintf("%.4f to %.4f", 
                        min(correlation_array, na.rm = TRUE), 
                        max(correlation_array, na.rm = TRUE)), "\n")
  
  return(invisible(list(
    correlation_array = correlation_array,
    mean_correlations = mean_correlations,
    sd_correlations = sd_correlations
  )))
}


run_analysis <- function() {
  base_input_dir <- "Results"
  
  tryCatch({
    results <- dimtest(base_input_dir)
    cat("Analysis completed successfully!\n")
    return(results)
  }, error = function(e) {
    cat("Analysis failed with error:", e$message, "\n")
    return(NULL)
  })
}


results <- run_analysis()