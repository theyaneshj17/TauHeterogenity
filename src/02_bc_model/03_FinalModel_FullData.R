# Load required libraries
library(coda)
library(MASS)
library(pROC)
library(corrplot)
library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
latent_dim <-5 #Optimal dimension selected per cluster

# Function to source all R files in a folder
sourceEntireFolder <- function(folderName, verbose=FALSE, showWarnings=TRUE) {
  files <- list.files(folderName, full.names=TRUE)
  files <- files[grepl("\\.[rR]$", files)]
  if (!length(files) && showWarnings)
    warning("No R files in ", folderName)
  
  for (f in files) {
    if (verbose) 
      cat("sourcing: ", f, "\n") 
    try(source(f, local=FALSE, echo=FALSE), silent=!verbose)
  }
  return(invisible(NULL))
}

sourceEntireFolder("code", verbose = FALSE, showWarnings = TRUE)


analyze_cluster_connectivity_bc <- function(matrix_list, cluster_assignments, output_base_dir, latent_dim, nscan=10000, burn=500) {

  output_dir <- paste0(output_base_dir, "/K=", latent_dim)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  

  log_file <- file(paste0(output_dir, "/log_latent_dim_", latent_dim, ".txt"), open="wt")
  sink(log_file, append=TRUE, type="output")
  sink(log_file, append=TRUE, type="message")
  

  cat("=== Running Analysis for Latent Dim =", latent_dim, " ===\n")
  cat("Start Time:", Sys.time(), "\n")
  cat("Output Directory:", output_dir, "\n")
  cat("nscan:", nscan, ", burn:", burn, "\n\n")
  
  # Split matrices by clusters
  cluster_data <- list()
  for(i in 0:2) {
    cluster_indices <- which(cluster_assignments == i)
    cat(sprintf("\nCluster %d: %d subjects\n", i, length(cluster_indices)))
    cluster_data[[i+1]] <- matrix_list[cluster_indices]
  }
  
  # Function to run BC model for each cluster
  run_bc_model <- function(cluster_matrices, cluster_num) {
    cat(sprintf("\n\nFitting BC model for Cluster %d", cluster_num))
    cat("\nNumber of subjects in cluster:", length(cluster_matrices))
    
    # Run BC model
    bc_result <- bc(
      X = cluster_matrices,
      W = NULL,
      K = latent_dim,  # Number of latent dimensions for connectivity modeling
      seed = 1,
      nscan = nscan,
      burn = burn,
      odens = 25,
      print = TRUE,
      gof = FALSE,
      plot = FALSE
    )
 
  saveRDS(bc_result, file = paste0(output_dir, "/bc_model_cluster_", cluster_num, ".rds"))
  cat(sprintf("\nBC model for Cluster %d saved to file: %s", cluster_num, paste0(output_dir, "/bc_model_cluster_", cluster_num, ".rds")))
  
    return(bc_result)
  }
  
  # Run BC model for each cluster
  bc_results <- list()
  for(i in 1:3) {
    bc_results[[i]] <- run_bc_model(cluster_data[[i]], i-1)
  }
  
  # Extract connectivity matrices from BC results
  connectivity_matrices <- list()
  for(i in 1:3) {
    connectivity_matrices[[i]] <- bc_results[[i]]$UVPM
  }
  
  # Generate heatmaps for each cluster's connectivity
  for(i in 1:3) {
    pdf(file = paste0(output_dir, "/cluster_", i-1, "_bc_connectivity.pdf"))
    corrplot(connectivity_matrices[[i]],
             method = "color",
             col = colorRampPalette(c("blue", "white", "red"))(200),
             tl.pos = "lt",
             tl.col = "black",
             tl.cex = 0.3,
             addgrid.col = "white",
             number.cex = 0.1,
             cl.pos = "b",
             is.corr = FALSE,
             main = paste("Cluster", i-1, "BC Model Estimated Connectivity"))
    dev.off()
  }
  
  # Compute and visualize differences between clusters
  diff_matrices <- list()
  diff_matrices[[1]] <- connectivity_matrices[[2]] - connectivity_matrices[[1]]  # Cluster 1 - Cluster 0
  diff_matrices[[2]] <- connectivity_matrices[[3]] - connectivity_matrices[[2]]  # Cluster 2 - Cluster 1
  diff_matrices[[3]] <- connectivity_matrices[[3]] - connectivity_matrices[[1]]  # Cluster 2 - Cluster 0
  
  # Generate difference heatmaps
  comparisons <- list(c(1,0), c(2,1), c(2,0))
  for(i in 1:3) {
    pdf(file = paste0(output_dir, "/cluster_", comparisons[[i]][1], "_vs_", 
                      comparisons[[i]][2], "_bc_diff.pdf"))
    corrplot(diff_matrices[[i]],
             method = "color",
             col = colorRampPalette(c("blue", "white", "red"))(200),
             tl.pos = "lt",
             tl.col = "black",
             tl.cex = 0.3,
             addgrid.col = "white",
             number.cex = 0.1,
             cl.pos = "b",
             is.corr = FALSE,
             main = paste("BC Model Connectivity Difference:\nCluster", 
                          comparisons[[i]][1], "-", comparisons[[i]][2]))
    dev.off()
  }
  

  summary_stats <- data.frame(
    Comparison = c("Cluster 1-0", "Cluster 2-1", "Cluster 2-0"),
    Mean_Diff = sapply(diff_matrices, mean),
    Max_Diff = sapply(diff_matrices, max),
    Min_Diff = sapply(diff_matrices, min)
  )
  

  write.csv(summary_stats, 
            file = paste0(output_dir, "/bc_cluster_comparison_summary.csv"), 
            row.names = FALSE)
  

  cat("End Time:", Sys.time(), "\n")
  

  sink()
  close(log_file)
  

  return(list(
    bc_results = bc_results,
    connectivity_matrices = connectivity_matrices,
    difference_matrices = diff_matrices,
    summary_stats = summary_stats
  ))
}


run_analysis_batch <- function(latent_dim, nscan, burn, output_base_dir) {
  # Load input data
  data <- readRDS("X.rds") # List of subject-level functional connectivity matrices
  pet_data <- read.csv("SubjectMapping.csv") # Subject-to-subtype cluster assignments from prior subtype discovery
  
  # Map subject IDs to cluster labels
  cluster_mapping <- setNames(pet_data$PredictedLabel3, pet_data$SubjectID)
  matched_labels <- cluster_mapping[names(data)]
  
  # Filter out NA values
  if (any(is.na(matched_labels))) {
    cat("Found", sum(is.na(matched_labels)), "NA values in matched_labels. Removing them.\n")
    valid_indices <- !is.na(matched_labels)
    data <- data[valid_indices]
    matched_labels <- matched_labels[valid_indices]
  }
  
  # Run the analysis for the current latent dimension
  output_dir <- paste0(output_base_dir, "/K=", latent_dim)
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  
  results <- analyze_cluster_connectivity_bc(
    data, 
    matched_labels, 
    output_dir, 
    latent_dim, 
    nscan, 
    burn
  )
}


nscan <- 50000           # Number of MCMC samples
burn <- 5000             # Number of burn-in samples
output_base_dir <- "Results"

run_analysis_batch(latent_dim, nscan, burn, output_base_dir)
