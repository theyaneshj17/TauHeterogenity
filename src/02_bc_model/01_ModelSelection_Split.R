# Load required libraries
library(coda)
library(MASS)
library(pROC)
library(corrplot)
library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
latent_dim <- as.numeric(args[1])

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
  
  
  
  run_bc_model <- function(cluster_matrices, cluster_num) {
    set.seed(47) 
    
    # Split into 80% train and 20% test
    total_subjects <- length(cluster_matrices)
    train_indices <- sample(seq_len(total_subjects), size = floor(0.8 * total_subjects))
    test_indices <- setdiff(seq_len(total_subjects), train_indices)
    
    train_cluster_matrices <- cluster_matrices[train_indices]
    test_cluster_matrices <- cluster_matrices[test_indices]
    
    # Save test data for later evaluation
    saveRDS(test_cluster_matrices, file = paste0(output_dir, "/test_data_cluster_", cluster_num, ".rds"))
    cat(sprintf("\nSaved test data for Cluster %d to file: %s", cluster_num, paste0(output_dir, "/test_data_cluster_", cluster_num, ".rds")))
    
    # Run BC model on training data
    cat(sprintf("\n\nFitting BC model for Cluster %d (Training Data)", cluster_num))
    cat("\nNumber of subjects in training cluster:", length(train_cluster_matrices))
    
    bc_result <- bc(
      X = train_cluster_matrices,
      W = NULL,
      K = latent_dim,
      seed = 1,
      nscan = nscan,
      burn = burn,
      odens = 25,
      print = TRUE,
      gof = FALSE,
      plot = FALSE
    )
    
    # Save the BC model result
    saveRDS(bc_result, file = paste0(output_dir, "/bc_model_cluster_", cluster_num, ".rds"))
    cat(sprintf("\nBC model for Cluster %d saved to file: %s", cluster_num, paste0(output_dir, "/bc_model_cluster_", cluster_num, ".rds")))
    
    return(list(
      bc_result = bc_result,
      test_data = test_cluster_matrices
    ))
  }
  
  
  # Run BC model for each cluster
  bc_results <- list()
  for(i in 1:3) {
    run_bc_model(cluster_data[[i]], i-1)
  }
  
  
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
