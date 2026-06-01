#Perform for both 9 and 17 groups in parallel sessions
#set working directory for directory with EINS RDS file
setwd()
#set seed
set.seed(47)
library(randomForest)
library(caret)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tibble)
library(gridExtra)
library(xgboost)
library(e1071)
library(glmnet)

#open RDS file
EINS_CCLE <- readRDS(file = "CCLE_Super_Eval_17.rds")

##########################################################
# HELPER FUNCTIONS
##########################################################

calculate_pca <- function(EINS_obj, omics_name, sample_names, n_comp = 9) {
  data_mat <- EINS_obj$Preprocessed_Omics$Matched_Data[[omics_name]]
  if (is.null(data_mat)) return(NULL)
  if (nrow(data_mat) != length(sample_names)) data_mat <- t(data_mat)
  
  col_var <- apply(data_mat, 2, var, na.rm = TRUE)
  data_mat <- data_mat[, col_var > 0, drop = FALSE]
  
  pca_res <- prcomp(data_mat, center = TRUE, scale. = TRUE, rank. = n_comp)
  pca_scores <- pca_res$x[match(sample_names, rownames(pca_res$x)), , drop = FALSE]
  
  return(pca_scores)
}

extract_embedding <- function(EINS_obj, method, sample_names, k = 9) {
  coord_path <- list(
    MoCluster = c("MoCluster", paste0("Factors_", k)),
    LRAcluster = c("LRAcluster", paste0("Factors_", k)),
    GAUDI = c("GAUDI", "Kmeans", paste0("Factors_", k)),
    COCA = c("COCA", paste0("Factors_", k)),
    SNF = c("SNF", "euclidean squared", paste0("Factors_", k)),
    MCIA = c("MCIA", paste0("Factors_", k)),
    MOFA = c("MOFA", paste0("Factors_", k))
  )
  
  # Try 1: Extract from CoordData if available
  if (method %in% names(coord_path)) {
    emb <- EINS_obj$Multi_Omics$CoordData
    for (key in coord_path[[method]]) {
      if (!is.null(emb)) emb <- emb[[key]]
    }
    if (!is.null(emb)) {
      emb <- emb[match(sample_names, rownames(emb)), , drop = FALSE]
      if (!any(is.na(emb))) {
        cat("  ", method, "- using CoordData embeddings\n")
        return(emb)
      }
    }
  }
  
  # Try 2: For NMF methods, use W matrix (sample factor loadings)
  if (method %in% c("iNMF", "jNMF")) {
    fit_res <- EINS_obj$Multi_Omics$Fit[[method]][[paste0("Clusters_", k)]]
    
    if (!is.null(fit_res) && !is.null(fit_res$W)) {
      W_matrix <- fit_res$W
      
      if (!is.null(dim(W_matrix))) {
        W_matrix <- W_matrix[match(sample_names, rownames(W_matrix)), , drop = FALSE]
        
        if (!any(is.na(W_matrix))) {
          colnames(W_matrix) <- paste0(method, "_Factor_", 1:ncol(W_matrix))
          cat("  ", method, "- using W matrix from Fit results\n")
          return(W_matrix)
        }
      }
    }
  }
  
  # Try 3: Use cluster assignments as one-hot encoding (fallback for methods with only clustering)
  cluster_path <- list(
    MoCluster = c("MoCluster", paste0("Clusters_", k)),
    LRAcluster = c("LRAcluster", paste0("Clusters_", k)),
    GAUDI = c("GAUDI", "Kmeans", paste0("Clusters_", k)),
    COCA = c("COCA", paste0("Clusters_", k)),
    SNF = c("SNF", "euclidean squared", paste0("Clusters_", k)),
    iNMF = c("iNMF", paste0("Clusters_", k)),
    jNMF = c("jNMF", paste0("Clusters_", k))
  )
  
  if (method %in% names(cluster_path)) {
    clust_res <- EINS_obj$Multi_Omics$ClusterRes
    for (key in cluster_path[[method]]) {
      if (!is.null(clust_res)) clust_res <- clust_res[[key]]
    }
    
    if (!is.null(clust_res)) {
      clust_res <- clust_res[match(sample_names, rownames(clust_res)), , drop = FALSE]
      cluster_labels <- as.numeric(clust_res$Cluster)
      
      n_samples <- length(cluster_labels)
      n_clusters <- k
      one_hot <- matrix(0, nrow = n_samples, ncol = n_clusters)
      
      for (i in 1:n_samples) {
        if (!is.na(cluster_labels[i]) && cluster_labels[i] <= n_clusters) {
          one_hot[i, cluster_labels[i]] <- 1
        }
      }
      
      rownames(one_hot) <- sample_names
      colnames(one_hot) <- paste0(method, "_Cluster_", 1:n_clusters)
      cat("  ", method, "- using one-hot encoded clusters (fallback)\n")
      return(one_hot)
    }
  }
  
  cat("  ", method, "- NOT FOUND\n")
  return(NULL)
}

calculate_ensemble_mds <- function(EINS_obj, sample_names, k = 9) {
  ens_clust <- EINS_obj$Ensemble$Samples$MDS$ClusterRes[[paste0("Clusters_", k)]]
  if (is.null(ens_clust)) return(NULL)
  
  ens_clust <- ens_clust[match(sample_names, rownames(ens_clust)), , drop = FALSE]
  labels <- as.numeric(ens_clust$Cluster)
  
  co_mat <- outer(labels, labels, "==") * 1
  emb <- cmdscale(as.dist(1 - co_mat), k = k)
  rownames(emb) <- sample_names
  
  return(emb)
}

calculate_ensemble_cca <- function(EINS_obj, sample_names, k = 9, reference_method = "MoCluster") {
  mcia <- EINS_obj$Multi_Omics$CoordData$MCIA[[paste0("Factors_", k)]]
  mofa <- EINS_obj$Multi_Omics$CoordData$MOFA[[paste0("Factors_", k)]]
  mocluster <- EINS_obj$Multi_Omics$CoordData$MoCluster[[paste0("Factors_", k)]]
  lracluster <- EINS_obj$Multi_Omics$CoordData$LRAcluster[[paste0("Factors_", k)]]
  gaudi <- EINS_obj$Multi_Omics$CoordData$GAUDI$Kmeans[[paste0("Factors_", k)]]
  coca <- EINS_obj$Multi_Omics$CoordData$COCA[[paste0("Clusters_", k)]]
  snf <- EINS_obj$Multi_Omics$CoordData$SNF$`euclidean squared`[[paste0("Clusters_", k)]]
  embeddings <- list()
  
  if (!is.null(mcia)) {
    mcia <- mcia[match(sample_names, rownames(mcia)), , drop = FALSE]
    if (!any(is.na(mcia))) embeddings[["MCIA"]] <- as.matrix(mcia)
  }
  
  if (!is.null(mofa)) {
    mofa <- mofa[match(sample_names, rownames(mofa)), , drop = FALSE]
    if (!any(is.na(mofa))) embeddings[["MOFA"]] <- as.matrix(mofa)
  }
  
  if (!is.null(mocluster)) {
    mocluster <- mocluster[match(sample_names, rownames(mocluster)), , drop = FALSE]
    if (!any(is.na(mocluster))) embeddings[["MoCluster"]] <- as.matrix(mocluster)
  }
  
  if (!is.null(lracluster)) {
    lracluster <- lracluster[match(sample_names, rownames(lracluster)), , drop = FALSE]
    if (!any(is.na(lracluster))) embeddings[["LRAcluster"]] <- as.matrix(lracluster)
  }
  
  if (!is.null(gaudi)) {
    gaudi <- gaudi[match(sample_names, rownames(gaudi)), , drop = FALSE]
    if (!any(is.na(gaudi))) embeddings[["GAUDI"]] <- as.matrix(gaudi)
  }
  
  if (!is.null(coca)) {
    coca <- coca[match(sample_names, rownames(coca)), , drop = FALSE]
    if (!any(is.na(coca))) embeddings[["COCA"]] <- as.matrix(coca)
  }
  
  if (!is.null(snf)) {
    snf <- snf[match(sample_names, rownames(snf)), , drop = FALSE]
    if (!any(is.na(snf))) embeddings[["SNF"]] <- as.matrix(snf)
  }
  
  iNMF_feat <- extract_embedding(EINS_obj, "iNMF", sample_names, k)
  if (!is.null(iNMF_feat)) embeddings[["iNMF"]] <- iNMF_feat
  
  jNMF_feat <- extract_embedding(EINS_obj, "jNMF", sample_names, k)
  if (!is.null(jNMF_feat)) embeddings[["jNMF"]] <- jNMF_feat
  
  if (length(embeddings) < 2) return(NULL)
  
  if (reference_method %in% names(embeddings)) {
    ref_idx <- which(names(embeddings) == reference_method)
    embeddings <- c(embeddings[ref_idx], embeddings[-ref_idx])
  }
  
  reference <- as.matrix(embeddings[[1]])
  n_samples <- nrow(reference)
  
  aligned_list <- list()
  aligned_list[[1]] <- reference
  
  for (i in 2:length(embeddings)) {
    current_emb <- as.matrix(embeddings[[i]])
    
    if (nrow(current_emb) != n_samples) {
      aligned_list[[i]] <- current_emb
      next
    }
    
    tryCatch({
      n_cca_comp <- min(ncol(reference), ncol(current_emb), n_samples - 1)
      
      cca_result <- cancor(x = as.matrix(reference[, 1:n_cca_comp, drop = FALSE]),
                           y = as.matrix(current_emb[, 1:n_cca_comp, drop = FALSE]))
      
      if (!is.null(cca_result$ycoef) && is.matrix(cca_result$ycoef)) {
        n_canon <- min(ncol(cca_result$ycoef), n_cca_comp)
        aligned_current <- as.matrix(current_emb[, 1:n_cca_comp, drop = FALSE]) %*%
          as.matrix(cca_result$ycoef[, 1:n_canon, drop = FALSE])
        aligned_list[[i]] <- aligned_current
      } else {
        aligned_list[[i]] <- current_emb
      }
    }, error = function(e) {
      aligned_list[[i]] <<- current_emb
    })
  }
  
  aligned_embeddings <- do.call(cbind, aligned_list)
  rownames(aligned_embeddings) <- sample_names
  colnames(aligned_embeddings) <- paste0("CCA_", 1:ncol(aligned_embeddings))
  
  return(aligned_embeddings)
}

cross_validation <- function(data_mat, labels, method_name, 
                             classifiers = c("RandomForest", "SVM"),
                             n_repeats = 100, fold_size = 0.2, seed = 42) {
  set.seed(seed = seed)
  
  n_samples <- nrow(data_mat)
  test_size <- floor(n_samples * fold_size)
  
  # Pre-generate all train/test splits so all classifiers use identical data
  all_splits <- vector("list", n_repeats)
  for (i in 1:n_repeats) {
    test_idx <- sample(1:n_samples, size = test_size, replace = FALSE)
    train_idx <- setdiff(1:n_samples, test_idx)
    all_splits[[i]] <- list(train = train_idx, test = test_idx)
  }
  
  results_list <- list()
  
  for (clf in classifiers) {
    all_acc <- numeric(n_repeats)
    all_f1 <- numeric(n_repeats)
    for (i in 1:n_repeats) {
      if(i %% 10 == 0) {
        cat("CV iteration:", i, "for classifier:", clf, "\n")
      }
      
      train_idx <- all_splits[[i]]$train
      test_idx <- all_splits[[i]]$test
      
      train_x <- data_mat[train_idx, , drop = FALSE]
      train_y <- labels[train_idx]
      test_x <- data_mat[test_idx, , drop = FALSE]
      test_y <- labels[test_idx]
      
      tryCatch({
        if (clf == "RandomForest") {
          model <- randomForest(
            x = train_x,
            y = train_y,
            ntree = 1000,
            importance = FALSE
          )
          pred <- stats::predict(object = model, newdata = test_x)
          
        } else if (clf == "XGBoost") {
          label_map <- setNames(0:(nlevels(train_y)-1), levels(train_y))
          train_labels_num <- label_map[as.character(train_y)]
          
          train_x_mat <- as.matrix(train_x)
          test_x_mat <- as.matrix(test_x)
          
          dtrain <- xgb.DMatrix(data = train_x_mat, label = train_labels_num)
          dtest <- xgb.DMatrix(data = test_x_mat)
          
          params <- list(
            objective = "multi:softmax",
            num_class = nlevels(train_y),
            eta = 0.3,
            max_depth = 6,
            subsample = 0.8
          )
          
          model <- xgb.train(
            params = params,
            data = dtrain,
            nrounds = 100,
            verbose = 0
          )
          
          pred_num <- stats::predict(model, dtest)
          pred <- factor(names(label_map)[match(pred_num, label_map)], levels = levels(train_y))
          
        } else if (clf == "SVM") {
          gamma_val <- 1 / ncol(train_x)
          model <- svm(
            x = train_x,
            y = train_y,
            kernel = "radial",
            cost = 1,
            gamma = gamma_val
          )
          pred <- stats::predict(model, newdata = test_x)
          
        } else if (clf == "LogisticRegression") {
          train_x_mat <- as.matrix(train_x)
          test_x_mat <- as.matrix(test_x)
          
          cv_fit <- cv.glmnet(
            x = train_x_mat,
            y = train_y,
            family = "multinomial",
            alpha = 0,
            nfolds = 5
          )
          pred <- stats::predict(cv_fit, newx = test_x_mat, s = "lambda.min", type = "class")
          pred <- factor(pred[,1], levels = levels(train_y))
          
        } else if (clf == "LASSO") {
          train_x_mat <- as.matrix(train_x)
          test_x_mat <- as.matrix(test_x)
          
          cv_fit <- cv.glmnet(
            x = train_x_mat,
            y = train_y,
            family = "multinomial",
            alpha = 1,
            nfolds = 5
          )
          pred <- stats::predict(cv_fit, newx = test_x_mat, s = "lambda.min", type = "class")
          pred <- factor(pred[,1], levels = levels(train_y))
          
        } else if (clf == "LASSO08") {
          train_x_mat <- as.matrix(train_x)
          test_x_mat <- as.matrix(test_x)
          
          cv_fit <- cv.glmnet(
            x = train_x_mat,
            y = train_y,
            family = "multinomial",
            alpha = 0.8,
            nfolds = 5
          )
          pred <- stats::predict(cv_fit, newx = test_x_mat, s = "lambda.min", type = "class")
          pred <- factor(pred[,1], levels = levels(train_y))
        }
        
        cm <- confusionMatrix(pred, test_y)
        all_acc[i] <- cm$overall['Accuracy']
        all_f1[i] <- mean(cm$byClass[, 'F1'], na.rm = TRUE)
        
      }, error = function(e) {
        all_acc[i] <<- NA
        all_f1[i] <<- NA
      })
    }
    
    all_acc <- all_acc[!is.na(all_acc)]
    all_f1 <- all_f1[!is.na(all_f1)]
    
    results_list[[clf]] <- list(
      method = method_name,
      classifier = clf,
      accuracy = mean(all_acc),
      acc_sd = sd(all_acc),
      f1_macro = mean(all_f1),
      f1_sd = sd(all_f1),
      all_accuracies = all_acc,
      all_f1_scores = all_f1
    )
  }
  
  return(results_list)
}

plot_f1_methods <- function(data = combined_results,
                            classifier = "RandomForest",
                            methods = NULL,
                            ensemble_only = "Ensemble_CCA",
                            save_plot = FALSE,
                            filename = NULL) {
  
  valid_classifiers <- c("RandomForest", "XGBoost", "SVM", "LogisticRegression", "LASSO", "LASSO08")
  if (!classifier %in% valid_classifiers) {
    stop("Invalid classifier. Choose from: ", paste(valid_classifiers, collapse = ", "))
  }
  
  filtered_data <- data %>%
    filter(Classifier == classifier)
  
  if (!is.null(ensemble_only)) {
    filtered_data <- filtered_data %>%
      filter(!grepl("Ensemble", Method) | Method %in% ensemble_only)
  }
  
  if (!is.null(methods)) {
    filtered_data <- filtered_data %>%
      filter(Method %in% methods)
    
    if (nrow(filtered_data) == 0) {
      stop("No data found for specified methods: ", paste(methods, collapse = ", "))
    }
  }
  
  plot_data <- filtered_data %>%
    group_by(Method) %>%
    summarise(
      Mean_F1 = mean(F1_Macro),
      SD_F1 = mean(F1_SD),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_F1)) %>%
    mutate(Method = factor(Method, levels = Method))
  
  subtitle_parts <- paste0(classifier, " classifier")
  if (!is.null(ensemble_only)) {
    subtitle_parts <- paste0(subtitle_parts, " | Ensemble: ", paste(ensemble_only, collapse = ", "))
  } else {
    subtitle_parts <- paste0(subtitle_parts, " | All ensemble methods included")
  }
  if (!is.null(methods)) {
    subtitle_parts <- paste0(subtitle_parts, " | ", length(methods), " methods selected")
  }
  subtitle_parts <- paste0(subtitle_parts, " | Error bars: F1 SD")
  
  p <- ggplot(plot_data, aes(x = Method, y = Mean_F1)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_errorbar(aes(ymin = Mean_F1 - SD_F1, ymax = Mean_F1 + SD_F1),
                  width = 0.3, linewidth = 0.8, color = "darkred") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "F1 Macro Performance by Method (Sorted Highest to Lowest)",
      subtitle = subtitle_parts,
      x = "Method",
      y = "F1 Macro Score"
    ) +
    ylim(0, max(plot_data$Mean_F1 + plot_data$SD_F1) * 1.05)
  
  if (save_plot) {
    if (is.null(filename)) {
      filename <- paste0("f1_methods_", classifier, ".pdf")
    }
    ggsave(filename, p, width = 12, height = 8)
    cat("Plot saved to:", filename, "\n")
  }
  
  return(p)
}

plot_f1_classifiers <- function(data = combined_results,
                                method = "Ensemble_CCA",
                                save_plot = FALSE,
                                filename = NULL) {
  
  if (!method %in% unique(data$Method)) {
    stop("Method '", method, "' not found in data. Available methods:\n",
         paste(unique(data$Method), collapse = ", "))
  }
  
  plot_data <- data %>%
    filter(Method == method) %>%
    group_by(Classifier) %>%
    summarise(
      Mean_F1 = mean(F1_Macro),
      SD_F1 = mean(F1_SD),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_F1)) %>%
    mutate(Classifier = factor(Classifier, levels = Classifier))
  
  p <- ggplot(plot_data, aes(x = Classifier, y = Mean_F1)) +
    geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
    geom_errorbar(aes(ymin = Mean_F1 - SD_F1, ymax = Mean_F1 + SD_F1),
                  width = 0.3, linewidth = 0.8, color = "darkred") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste0("F1 Macro Performance: ", method),
      subtitle = "Comparison across all classifiers | Sorted highest to lowest | Error bars: F1 SD",
      x = "Classifier",
      y = "F1 Macro Score"
    ) +
    ylim(0, max(plot_data$Mean_F1 + plot_data$SD_F1) * 1.05)
  
  if (save_plot) {
    if (is.null(filename)) {
      filename <- paste0("f1_classifiers_", gsub(" ", "_", method), ".pdf")
    }
    ggsave(filename, p, width = 10, height = 7)
    cat("Plot saved to:", filename, "\n")
  }
  
  return(p)
}

plot_f1_all_methods_classifiers <- function(data = combined_results,
                                            save_plot = FALSE,
                                            filename = NULL) {
  
  single_omics <- c("Methylation", "miRNA", "Transcriptomics", "Proteomics", "Metabolomics")
  
  plot_data <- data %>%
    mutate(
      Category = case_when(
        Method %in% single_omics ~ "Single-Omics",
        grepl("Ensemble", Method) ~ "Ensemble",
        TRUE ~ "Multi-Omics"
      ),
      Method_Classifier = paste0(Method, " (", Classifier, ")")
    ) %>%
    arrange(desc(F1_Macro))
  
  plot_data$Method_Classifier <- factor(plot_data$Method_Classifier, 
                                        levels = unique(plot_data$Method_Classifier))
  
  category_colors <- c(
    "Single-Omics" = "#E69F00",
    "Multi-Omics" = "#56B4E9", 
    "Ensemble" = "#009E73"
  )
  
  p <- ggplot(plot_data, aes(x = Method_Classifier, y = F1_Macro, fill = Category)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_errorbar(aes(ymin = F1_Macro - F1_SD, ymax = F1_Macro + F1_SD),
                  width = 0.4, linewidth = 0.6, color = "black", alpha = 0.6) +
    scale_fill_manual(values = category_colors) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10)
    ) +
    labs(
      title = "F1 Macro Performance: All Methods × All Classifiers",
      subtitle = "Sorted by performance | Error bars: F1 SD | Colors indicate method category",
      x = "Method (Classifier)",
      y = "F1 Macro Score",
      fill = "Method Category"
    ) +
    ylim(0, max(plot_data$F1_Macro + plot_data$F1_SD) * 1.05)
  
  if (save_plot) {
    if (is.null(filename)) {
      filename <- "f1_all_methods_classifiers.pdf"
    }
    ggsave(filename, p, width = 16, height = 8)
    cat("Plot saved to:", filename, "\n")
  }
  
  return(p)
}

plot_f1_methods_averaged <- function(data = combined_results,
                                     save_plot = FALSE,
                                     filename = NULL) {
  
  single_omics <- c("Methylation", "miRNA", "Transcriptomics", "Proteomics", "Metabolomics")
  
  plot_data <- data %>%
    mutate(
      Category = case_when(
        Method %in% single_omics ~ "Single-Omics",
        grepl("Ensemble", Method) ~ "Ensemble",
        TRUE ~ "Multi-Omics"
      )
    ) %>%
    group_by(Method, Category) %>%
    summarise(
      Mean_F1 = mean(F1_Macro),
      SD_F1 = sd(F1_Macro),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_F1)) %>%
    mutate(Method = factor(Method, levels = Method))
  
  category_colors <- c(
    "Single-Omics" = "#E69F00",
    "Multi-Omics" = "#56B4E9", 
    "Ensemble" = "#009E73"
  )
  print(plot_data, n = 24)
  
  p <- ggplot(plot_data, aes(x = Method, y = Mean_F1, fill = Category)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_errorbar(aes(ymin = Mean_F1 - SD_F1, ymax = Mean_F1 + SD_F1),
                  width = 0.4, linewidth = 0.6, color = "black", alpha = 0.6) +
    scale_fill_manual(values = category_colors) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 16),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 16)
    ) +
    labs(
      title = "F1 Macro Performance by Method (Averaged Across Classifiers)",
      subtitle = "Sorted by performance | Error bars: SD across classifiers | Colors indicate method category",
      x = "Method",
      y = "F1 Macro Score (Mean)",
      fill = "Method Category"
    ) +
    ylim(0, max(plot_data$Mean_F1 + plot_data$SD_F1, na.rm = TRUE) * 1.05)
  
  if (save_plot) {
    if (is.null(filename)) {
      filename <- "f1_methods_averaged.pdf"
    }
    ggsave(filename, p, width = 12, height = 8)
    cat("Plot saved to:", filename, "\n")
  }
  
  return(p)
}

##########################################################
# MAIN ANALYSIS
##########################################################

ds_name <- "CCLE_17"
EINS_obj <- EINS_CCLE

##

true_labels <- as.factor(EINS_obj$Omics$Metadata$Proteomics$Site_Primary)
sample_names <- rownames(EINS_obj$Omics$Metadata$Proteomics)

n_classes <- nlevels(true_labels)
cat("Samples:", length(true_labels), "| Classes:", n_classes, "\n")

cat("Extracting representations...\n")
all_data <- list()

single_omics <- c("Methylation", "miRNA", "Transcriptomics", "Proteomics", "Metabolomics")
for (omics in single_omics) {
  all_data[[omics]] <- calculate_pca(EINS_obj, omics, sample_names, 9)
}

multiomics_methods <- c("MCIA","MOFA","MoCluster", "LRAcluster", "GAUDI", "COCA", "SNF", "iNMF", "jNMF")

for (method in multiomics_methods) {
  all_data[[method]] <- extract_embedding(EINS_obj, method, sample_names, 9)
}

all_data[["Ensemble_MDS"]] <- calculate_ensemble_mds(EINS_obj, sample_names, 9)
all_data[["Ensemble_CCA_MCIA"]] <- calculate_ensemble_cca(EINS_obj, sample_names, 9, "MCIA")
all_data[["Ensemble_CCA_MOFA"]] <- calculate_ensemble_cca(EINS_obj, sample_names, 9, "MOFA")
all_data[["Ensemble_CCA_MoCluster"]] <- calculate_ensemble_cca(EINS_obj, sample_names, 9, "MoCluster")
all_data[["Ensemble_CCA_LRAcluster"]] <- calculate_ensemble_cca(EINS_obj, sample_names, 9, "LRAcluster")
all_data[["Ensemble_CCA_GAUDI"]] <- calculate_ensemble_cca(EINS_obj, sample_names, 9, "GAUDI")
all_data[["Ensemble_CCA_COCA"]] <- calculate_ensemble_cca(EINS_obj, sample_names, 9, "COCA")
all_data[["Ensemble_CCA_SNF"]] <- calculate_ensemble_cca(EINS_obj, sample_names, 9, "SNF")
all_data[["Ensemble_CCA_iNMF"]] <- calculate_ensemble_cca(EINS_obj, sample_names, 9, "iNMF")
all_data[["Ensemble_CCA_jNMF"]] <- calculate_ensemble_cca(EINS_obj, sample_names, 9, "jNMF")

all_data <- all_data[!sapply(all_data, is.null)]

#Running cross-validation
results_list <- list()

for (method_name in names(all_data)) {
  data_mat <- all_data[[method_name]]
  
  if (any(is.na(data_mat))) {
    cat("  ", method_name, "- SKIPPED (contains NAs)\n")
    next
  }
  
  cat("  ", method_name, "...\n")
  
  res_classifiers <- 
    cross_validation(data_mat, true_labels, method_name, n_repeats = 20, fold_size = 0.2, seed = 21,
                     classifiers = c("RandomForest", "XGBoost", "SVM", "LogisticRegression", "LASSO", "LASSO08"))
  
  for (clf in names(res_classifiers)) {
    res <- res_classifiers[[clf]]
    results_list[[paste0(method_name, "_", clf)]] <- res
    cat("    ", clf, ": Acc =", round(res$accuracy*100, 1), "% ±", round(res$acc_sd*100, 1), 
        "%, F1 =", round(res$f1_macro*100, 1), "% ±", round(res$f1_sd*100, 1), "%\n")
  }
}

results_df <- do.call(rbind, lapply(results_list, function(x) {
  data.frame(
    Dataset = ds_name,
    N_Classes = n_classes,
    Method = x$method,
    Classifier = x$classifier,
    Accuracy = x$accuracy,
    Accuracy_SD = x$acc_sd,
    F1_Macro = x$f1_macro,
    F1_SD = x$f1_sd,
    stringsAsFactors = FALSE
  )
}))

combined_results <- results_df %>%
  mutate(
    Harmonic_Mean = 2 * (Accuracy * F1_Macro) / (Accuracy + F1_Macro),
    Geometric_Mean = sqrt(Accuracy * F1_Macro),
    Combined_SD = (Accuracy_SD + F1_SD) / 2
  )

cat("\n")
cat("TOP METHODS BY HARMONIC MEAN:\n")
top_harmonic <- combined_results %>%
  arrange(desc(Harmonic_Mean)) %>%
  head(20)
print(top_harmonic)

cat("\nANALYSIS COMPLETE!\n")

##########################################################
# F1 PLOT: Methods sorted by performance
##########################################################

dev.off()
combined_results2=combined_results
combined_results2$Method<-gsub("_iNMF", " (iNMF)", combined_results2$Method)
combined_results2$Method<-gsub("_jNMF", " (jNMF)", combined_results2$Method)
combined_results2$Method<-gsub("_COCA", " (COCA)", combined_results2$Method)
combined_results2$Method<-gsub("_SNF", " (SNF)", combined_results2$Method)
combined_results2$Method<-gsub("_LASSO", " (LASSO)", combined_results2$Method)
combined_results2$Method<-gsub("_LASSO08", " (LASSO08)", combined_results2$Method)
combined_results2$Method<-gsub("_LogisticRegression", " (Logistic Regression)", combined_results2$Method)
combined_results2$Method<-gsub("_LRAcluster", " (LRAcluster)", combined_results2$Method)
combined_results2$Method<-gsub("_GAUDI", " (GAUDI)", combined_results2$Method)
combined_results2$Method<-gsub("_MoCluster", " (MoCluster)", combined_results2$Method)
combined_results2$Method<-gsub("_MCIA", " (MCIA)", combined_results2$Method)
combined_results2$Method<-gsub("_MOFA", " (MOFA)", combined_results2$Method)

p_f1_avg_no_can <- plot_f1_methods_averaged(combined_results2)

