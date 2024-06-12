#' Run pseudotime estimation and ML pipeline
#' 
#' @param seurat_path Path to the Seurat object file (.rds)
#' @return List with results including RMSE, R-squared, top genes, and fitted models
#'
run_pseudotime_pipeline <- function(seurat_path) {
  library(Seurat)
  library(monocle3)
  library(caret)
  library(xgboost)
  
  # Load Seurat object
  seu <- readRDS(seurat_path)
  
  # Transform into CDS for monocle 
  cds <- SeuratWrappers::as.cell_data_set(seu)
  cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(seu[["RNA"]])
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds <- order_cells(cds)
  
  # Plot pseudotime
  plot_cells(cds,
             color_cells_by = "pseudotime",
             label_groups_by_cluster = FALSE,
             label_leaves = FALSE,
             label_branch_points = FALSE, 
             cell_size = 0.9, 
             trajectory_graph_segment_size = 2)
  
  # Extract gene expression matrix and pseudotime values
  expr_matrix <- as.matrix(exprs(cds))
  pseudotime_values <- pseudotime(cds)
  
  # Combine into a data frame
  data <- data.frame(t(expr_matrix))
  data$pseudotime <- pseudotime_values
  
  # Split into training and testing sets
  set.seed(123)
  trainIndex <- createDataPartition(data$pseudotime, p = 0.8, list = FALSE)
  trainData <- data[trainIndex, ]
  testData <- data[-trainIndex, ]
  
  # Prepare data for xgboost
  x_train <- as.matrix(trainData[, -ncol(trainData)])
  y_train <- trainData$pseudotime
  x_test <- as.matrix(testData[, -ncol(testData)])
  y_test <- testData$pseudotime
  
  # Convert to xgboost DMatrix
  dtrain <- xgb.DMatrix(data = x_train, label = y_train)
  dtest <- xgb.DMatrix(data = x_test, label = y_test)
  
  # Train an xgboost model
  xgb_model <- xgboost(
    data = dtrain, 
    max_depth = 3, 
    nrounds = 100, 
    eta = 0.1, 
    objective = "reg:squarederror",
    verbose = 0
  )
  
  # Predict on the test set
  xgb_predictions <- predict(xgb_model, dtest)
  
  # Evaluate the model
  xgb_rmse <- sqrt(mean((y_test - xgb_predictions)^2))
  cat("xgboost RMSE:", xgb_rmse, "\n")
  
  # Extract feature importance scores
  importance <- xgb.importance(model = xgb_model)
  
  # Print the top 10 most important genes
  top_genes <- importance$Feature[1:10]
  cat("Top 10 important genes:\n")
  print(top_genes)
  
  # Calculate R-squared (coefficient of determination)
  SS_total <- sum((y_test - mean(y_test))^2)  # Total sum of squares
  SS_residual <- sum((y_test - xgb_predictions)^2)  # Residual sum of squares
  r2 <- 1 - (SS_residual / SS_total)
  cat("R-squared (coefficient of determination):", r2, "\n")
  
  # Fit models using specified genes
  
  cds_subset <- cds[rowData(cds)$gene_short_name %in% top_genes, ]
  gene_fits <- fit_models(cds_subset, model_formula_str = "~pseudotime")
  fit_coefs <- coefficient_table(gene_fits)
  time_terms <- fit_coefs %>% filter(term == "pseudotime")
  time_terms_filtered <- time_terms %>% filter(q_value < 0.05) %>%
    select(gene_short_name, term, q_value, estimate)
  
  # Return results as a list
  results <- list(
    xgb_rmse = xgb_rmse,
    r2 = r2,
    top_genes = top_genes,
    emb_time_terms = emb_time_terms_filtered,
    gene_fits = gene_fits
  )
  
  return(results)
}

# Example usage:
# results <- run_pseudotime_pipeline("/path/to/your/seurat_object.rds")
