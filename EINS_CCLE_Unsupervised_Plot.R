set.seed(47)
library(tidyverse)
library(aricode)
library(mclust)
library(ggplot2)
library(reshape2)
library(ggpubr)

#Changes filenames
EINS_CCLE3 <- readRDS(file = "CCLE_Unsuper_Eval_3.rds")
EINS_CCLE4 <- readRDS(file = "CCLE_Unsuper_Eval_4.rds")
EINS_CCLE5 <- readRDS(file = "CCLE_Unsuper_Eval_5.rds")
EINS_CCLE6 <- readRDS(file = "CCLE_Unsuper_Eval_6.rds")
EINS_CCLE7 <- readRDS(file = "CCLE_Unsuper_Eval_7.rds")
EINS_CCLE8 <- readRDS(file = "CCLE_Unsuper_Eval_8.rds")
EINS_CCLE9 <- readRDS(file = "CCLE_Unsuper_Eval_9.rds")

#########
#ARI/NMI analysis
#########

#########
#Functions
#########
evaluate_clustering_quality <- function(true_labels, predicted_labels) {
  ari <- adjustedRandIndex(true_labels, predicted_labels)
  nmi <- NMI(true_labels, predicted_labels)
  return(list(ARI = ari, NMI = nmi))
}

ARI_NMI_Clusters <- function(ClusterRes,
                             ClusterList,
                             Labels,
                             Samples){
  ARI_NMI_Results <- data.frame()
  for (num_groups in ClusterList) {
    
    true_labels <- Labels
    cluster_result <- ClusterRes[[num_groups]][order(match(rownames(ClusterRes[[num_groups]]), Samples)), , drop = FALSE]
    cluster_labels <- as.numeric(cluster_result$Cluster)
    clust_metrics <- evaluate_clustering_quality(true_labels, cluster_labels)
    
    ARI_NMI_Results <- rbind(ARI_NMI_Results, data.frame(
      Num_Groups = as.integer(strsplit(num_groups, "_")[[1]][2]),
      ARI_Ensemble = clust_metrics$ARI,
      NMI_Ensemble = clust_metrics$NMI
    ))
  }
  return(ARI_NMI_Results)
}

ARI_NMI_HClust <- function(Groups = 2:20,
                           Labels,
                           Samples,
                           Tree){
  ARI_NMI_Results <- data.frame()
  for (num_groups in Groups) {
    true_labels <- Labels
    cluster_labels <- cutree(Tree, k = num_groups)
    cluster_labels <- cluster_labels[order(match(names(cluster_labels), Samples))]
    clust_metrics <- evaluate_clustering_quality(true_labels, cluster_labels)
    
    ARI_NMI_Results <- rbind(ARI_NMI_Results, data.frame(
      Num_Groups = num_groups,
      ARI = clust_metrics$ARI,
      NMI = clust_metrics$NMI
    ))
  }
  return(ARI_NMI_Results)
}

results_combined <- function(results_list, labels = NULL, metrics = NULL,
                             facet_scales = "free_y", ncol = 2,
                             palette = "Dark2", line_sizes = 1) {
  if (!is.null(labels)) names(results_list) <- labels
  if (is.null(names(results_list))) names(results_list) <- paste0("Set ", seq_along(results_list))
  dataset_levels <- names(results_list)
  
  df_list <- Map(function(d, nm) { d$Dataset <- nm; d }, results_list, dataset_levels)
  all_df <- do.call(rbind, df_list)
  all_df$Dataset <- factor(all_df$Dataset, levels = dataset_levels)
  
  if (is.null(metrics)) metrics <- setdiff(names(all_df), c("Num_Groups", "Dataset"))
  long_df <- reshape2::melt(all_df, id.vars = c("Num_Groups", "Dataset"),
                            measure.vars = metrics, variable.name = "Metric", value.name = "Value")
  long_df$Num_Groups <- as.numeric(as.character(long_df$Num_Groups))
  long_df$Dataset <- factor(long_df$Dataset, levels = dataset_levels)
  return(long_df)
}

standardize_metrics <- function(d) {
  d <- as.data.frame(d, stringsAsFactors = FALSE, check.names = FALSE)
  nm <- names(d)
  nm <- sub("^ARI($|[_-].*)", "ARI", nm, ignore.case = TRUE)
  nm <- sub("^NMI($|[_-].*)", "NMI", nm, ignore.case = TRUE)
  names(d) <- nm
  d2 <- column_to_rownames(d, var = "Num_Groups")
  t <- t(d2)
  GM <- apply(t[, c(1:19)], 2, function(x) exp(mean(log(x))))
  d2$GeoMean <- GM
  d2$GeoMean[is.nan(d2$GeoMean)] <- 0
  d2 <- rownames_to_column(d2, var = "Num_Groups")
}

plot_results_combined <- function(Results,
                                  PointSize = 2,
                                  BaseSize = 14,
                                  Title,
                                  NumGroups,
                                  Lower = -0.02,
                                  Upper = 0.8
){
  plot <- ggplot(Results, aes(x = Num_Groups, y = Value, color = Dataset)) +
    geom_point(aes(color = Dataset, alpha = Dataset), size = PointSize) +
    geom_point(data = subset(Results, Dataset == "Ensemble"), aes(color = Dataset, alpha = Dataset), size = PointSize) +
    scale_alpha_manual(values = c(0.6, 1, 1, 1, 1, 1)) +
    theme_minimal(base_size = BaseSize) +
    labs(subtitle = Title) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.title.position = "plot") +
    scale_color_manual(values = c("#4d4d4d", "#0072b2", "#56b4e9", "#009e73", "#e69f00", "#f0e442")) +
    geom_line(aes(size = Dataset, alpha = Dataset)) +
    geom_line(data = subset(Results, Dataset == "Ensemble"), aes(size = Dataset, alpha = Dataset)) +
    scale_size_manual(values = c(2.8, 1, 1, 1, 1, 1)) +
    geom_vline(xintercept = NumGroups, color = "red", size = 1) +
    ylim(Lower, Upper)
  return(plot)
}

plot_results_methods_combined <- function(Results,
                                          PointSize = 2,
                                          BaseSize = 14,
                                          Title,
                                          NumGroups,
                                          Lower = -0.02,
                                          Upper = 0.8
){
  plot <- ggplot(Results, aes(x = Num_Groups, y = Value, color = Dataset)) +
    geom_point(aes(color = Dataset, alpha = Dataset), size = PointSize) +
    geom_point(data = subset(Results, Dataset == "Ensemble"), aes(color = Dataset, alpha = Dataset), size = PointSize) +
    scale_alpha_manual(values = c(0.6, 1, 1, 1, 1, 1, 1, 1)) +
    theme_minimal(base_size = BaseSize) +
    labs(subtitle = Title) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.title.position = "plot") +
    scale_color_manual(values = c("#4d4d4d", "#cc79a7", "#d55e00", "#999999", "#a6761d", "#1b9e77", "#7570b3", "#e7298a")) +
    geom_line(aes(size = Dataset, alpha = Dataset)) +
    geom_line(data = subset(Results, Dataset == "Ensemble"), aes(size = Dataset, alpha = Dataset)) +
    scale_size_manual(values = c(2.8, 1, 1, 1, 1, 1, 1, 1)) +
    geom_vline(xintercept = NumGroups, color = "red", size = 1) +
    ylim(Lower, Upper)
  return(plot)
}

###########
#EINS_CCLE3
###########
#True labels
CCLE_3_True <- as.factor(EINS_CCLE3$Omics$Metadata$Proteomics$Site_Primary)
levels(CCLE_3_True) <- c(1:3)
CCLE_3_True <- as.numeric(CCLE_3_True)

#True order
CCLE_3_Sam <- rownames(EINS_CCLE3$Omics$Metadata$Proteomics)

#Single-omics data prep
CCLE_3_Single <- list()
CCLE_3_Ens <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE3$Ensemble$Samples$CHC$ClusterRes,
                               ClusterList = names(EINS_CCLE3$Ensemble$Samples$CHC$ClusterRes),
                               Labels = CCLE_3_True,
                               Samples = CCLE_3_Sam)
CCLE_3_Ens <- standardize_metrics(CCLE_3_Ens)
CCLE_3_Single$Ensemble <- CCLE_3_Ens

CCLE_3_Omics<- names(EINS_CCLE3$Single_Omics$HClustTree)
for(Omics in CCLE_3_Omics){
  Res <- ARI_NMI_HClust(Tree = EINS_CCLE3$Single_Omics$HClustTree[[Omics]],
                        Labels = CCLE_3_True,
                        Samples = CCLE_3_Sam)
  Res <- standardize_metrics(Res)
  CCLE_3_Single[[Omics]] <- Res
}

#Multi-omics data prep
CCLE_3_Multi <- list()
CCLE_3_Multi$Ensemble <- CCLE_3_Ens

CCLE_3_Methods <- names(EINS_CCLE3$Multi_Omics$ClusterRes)
for(Method in CCLE_3_Methods){
  if(Method == "SNF"){
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE3$Multi_Omics$ClusterRes$SNF$`euclidean squared`,
                            ClusterList = names(EINS_CCLE3$Multi_Omics$ClusterRes$SNF$`euclidean squared`),
                            Labels = CCLE_3_True,
                            Samples = CCLE_3_Sam)
    Res <- standardize_metrics(Res)
    CCLE_3_Multi[[Method]] <- Res
  }else if(Method == "GAUDI"){
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE3$Multi_Omics$ClusterRes$GAUDI$Kmeans,
                            ClusterList = names(EINS_CCLE3$Multi_Omics$ClusterRes$GAUDI$Kmeans),
                            Labels = CCLE_3_True,
                            Samples = CCLE_3_Sam)
    Res <- standardize_metrics(Res)
    CCLE_3_Multi[[Method]] <- Res
  }else {
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE3$Multi_Omics$ClusterRes[[Method]],
                            ClusterList = names(EINS_CCLE3$Multi_Omics$ClusterRes[[Method]]),
                            Labels = CCLE_3_True,
                            Samples = CCLE_3_Sam)
    Res <- standardize_metrics(Res)
    CCLE_3_Multi[[Method]] <- Res
  }
}
CCLE_3_Single <- CCLE_3_Single[c("Ensemble", "Methylation", "miRNA", "Transcriptomics", "Proteomics", "Metabolomics")]
CCLE_3_Multi <- CCLE_3_Multi[c("Ensemble", "MoCluster", "LRAcluster", "GAUDI", "COCA", "SNF", "MCIA", "MOFA")]

#combine results in single dataframe
CCLE_3_Single_Res <- results_combined(
  results_list = CCLE_3_Single,
  metrics = c("ARI", "NMI", "GeoMean"),
  line_sizes = c(2.8, 1, 1, 1, 1, 1)
)

CCLE_3_Single_Res_ARI <- CCLE_3_Single_Res %>% filter(Metric == "ARI")
CCLE_3_Single_Res_NMI <- CCLE_3_Single_Res %>% filter(Metric == "NMI")
CCLE_3_Single_Res_Geo <- CCLE_3_Single_Res %>% filter(Metric == "GeoMean")

CCLE_3_Multi_Res <- results_combined(
  results_list = CCLE_3_Multi,
  metrics = c("ARI", "NMI", "GeoMean"),
  facet_scales = "free_y",
  ncol = 2,
  palette = "Dark2",
  line_sizes = c(2.8, 1, 1, 1, 1, 1,1,1)
)


CCLE_3_Multi_Res_ARI <- CCLE_3_Multi_Res %>% filter(Metric == "ARI")
CCLE_3_Multi_Res_NMI <- CCLE_3_Multi_Res %>% filter(Metric == "NMI")
CCLE_3_Multi_Res_Geo <- CCLE_3_Multi_Res %>% filter(Metric == "GeoMean")

#plot results
plot_CCLE_3_Single_ARI <- plot_results_combined(Results = CCLE_3_Single_Res_ARI, NumGroups = 3, Title = "ARI (# true groups = 3)", Lower = -0.05, Upper = 0.85)
plot_CCLE_3_Single_NMI <- plot_results_combined(Results = CCLE_3_Single_Res_NMI, NumGroups = 3, Title = "NMI (# true groups = 3)", Lower = -0.05, Upper = 0.85)
plot_CCLE_3_Single_Geo <- plot_results_combined(Results = CCLE_3_Single_Res_Geo, NumGroups = 3, Title = "Single-omics methods", Lower = -0.05, Upper = 0.8)
plot_CCLE_3_Single_Geo <- plot_CCLE_3_Single_Geo + theme(axis.text.x = element_blank())

plot_CCLE_3_Multi_ARI <- plot_results_methods_combined(Results = CCLE_3_Multi_Res_ARI, NumGroups = 3, Title = "ARI (# true groups = 3)", Lower = -0.05, Upper = 0.85)
plot_CCLE_3_Multi_NMI <- plot_results_methods_combined(Results = CCLE_3_Multi_Res_NMI, NumGroups = 3, Title = "NMI (# true groups = 3)", Lower = -0.05, Upper = 0.85)
plot_CCLE_3_Multi_Geo <- plot_results_methods_combined(Results = CCLE_3_Multi_Res_Geo, NumGroups = 3, Title = "Multi-omics methods", Lower = -0.05, Upper = 0.8)
plot_CCLE_3_Multi_Geo <- plot_CCLE_3_Multi_Geo + theme(axis.text.x = element_blank())

ARI_NMI_3 <- ggarrange(plot_CCLE_3_Single_ARI, plot_CCLE_3_Single_NMI, plot_CCLE_3_Multi_ARI, plot_CCLE_3_Multi_NMI, ncol = 4, nrow = 1)
Geo_3 <- ggarrange(plot_CCLE_3_Single_Geo, plot_CCLE_3_Multi_Geo, ncol = 2, nrow = 1)

###########
#EINS_CCLE4
###########
#True labels
CCLE_4_True <- as.factor(EINS_CCLE4$Omics$Metadata$Proteomics$Site_Primary)
levels(CCLE_4_True) <- c(1:4)
CCLE_4_True <- as.numeric(CCLE_4_True)

#True order
CCLE_4_Sam <- rownames(EINS_CCLE4$Omics$Metadata$Proteomics)

#Single-omics data prep
CCLE_4_Single <- list()
CCLE_4_Ens <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE4$Ensemble$Samples$CHC$ClusterRes,
                               ClusterList = names(EINS_CCLE4$Ensemble$Samples$CHC$ClusterRes),
                               Labels = CCLE_4_True,
                               Samples = CCLE_4_Sam)
CCLE_4_Ens <- standardize_metrics(CCLE_4_Ens)
CCLE_4_Single$Ensemble <- CCLE_4_Ens

CCLE_4_Omics<- names(EINS_CCLE4$Single_Omics$HClustTree)
for(Omics in CCLE_4_Omics){
  Res <- ARI_NMI_HClust(Tree = EINS_CCLE4$Single_Omics$HClustTree[[Omics]],
                        Labels = CCLE_4_True,
                        Samples = CCLE_4_Sam)
  Res <- standardize_metrics(Res)
  CCLE_4_Single[[Omics]] <- Res
}

#Multi-omics data prep
CCLE_4_Multi <- list()
CCLE_4_Multi$Ensemble <- CCLE_4_Ens

CCLE_4_Methods <- names(EINS_CCLE4$Multi_Omics$ClusterRes)
for(Method in CCLE_4_Methods){
  if(Method == "SNF"){
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE4$Multi_Omics$ClusterRes$SNF$`euclidean squared`,
                            ClusterList = names(EINS_CCLE4$Multi_Omics$ClusterRes$SNF$`euclidean squared`),
                            Labels = CCLE_4_True,
                            Samples = CCLE_4_Sam)
    Res <- standardize_metrics(Res)
    CCLE_4_Multi[[Method]] <- Res
  }else if(Method == "GAUDI"){
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE4$Multi_Omics$ClusterRes$GAUDI$Kmeans,
                            ClusterList = names(EINS_CCLE4$Multi_Omics$ClusterRes$GAUDI$Kmeans),
                            Labels = CCLE_4_True,
                            Samples = CCLE_4_Sam)
    Res <- standardize_metrics(Res)
    CCLE_4_Multi[[Method]] <- Res
  }else {
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE4$Multi_Omics$ClusterRes[[Method]],
                            ClusterList = names(EINS_CCLE4$Multi_Omics$ClusterRes[[Method]]),
                            Labels = CCLE_4_True,
                            Samples = CCLE_4_Sam)
    Res <- standardize_metrics(Res)
    CCLE_4_Multi[[Method]] <- Res
  }
}
CCLE_4_Single <- CCLE_4_Single[c("Ensemble", "Methylation", "miRNA", "Transcriptomics", "Proteomics", "Metabolomics")]
CCLE_4_Multi <- CCLE_4_Multi[c("Ensemble", "MoCluster", "LRAcluster", "GAUDI", "COCA", "SNF", "MCIA", "MOFA")]

#combine results in single dataframe
CCLE_4_Single_Res <- results_combined(
  results_list = CCLE_4_Single,
  metrics = c("ARI", "NMI", "GeoMean"),
  line_sizes = c(2.8, 1, 1, 1, 1, 1)
)

CCLE_4_Single_Res_ARI <- CCLE_4_Single_Res %>% filter(Metric == "ARI")
CCLE_4_Single_Res_NMI <- CCLE_4_Single_Res %>% filter(Metric == "NMI")
CCLE_4_Single_Res_Geo <- CCLE_4_Single_Res %>% filter(Metric == "GeoMean")

CCLE_4_Multi_Res <- results_combined(
  results_list = CCLE_4_Multi,
  metrics = c("ARI", "NMI", "GeoMean"),
  facet_scales = "free_y",
  ncol = 2,
  palette = "Dark2",
  line_sizes = c(2.8, 1, 1, 1, 1, 1,1,1)
)


CCLE_4_Multi_Res_ARI <- CCLE_4_Multi_Res %>% filter(Metric == "ARI")
CCLE_4_Multi_Res_NMI <- CCLE_4_Multi_Res %>% filter(Metric == "NMI")
CCLE_4_Multi_Res_Geo <- CCLE_4_Multi_Res %>% filter(Metric == "GeoMean")

#plot results
plot_CCLE_4_Single_ARI <- plot_results_combined(Results = CCLE_4_Single_Res_ARI, NumGroups = 4, Title = "ARI (# true groups = 4)", Lower = -0.05, Upper = 0.7)
plot_CCLE_4_Single_NMI <- plot_results_combined(Results = CCLE_4_Single_Res_NMI, NumGroups = 4, Title = "NMI (# true groups = 4)", Lower = -0.05, Upper = 0.7)
plot_CCLE_4_Single_Geo <- plot_results_combined(Results = CCLE_4_Single_Res_Geo, NumGroups = 4, Title = NULL, Lower = -0.05, Upper = 0.7)
plot_CCLE_4_Single_Geo <- plot_CCLE_4_Single_Geo + theme(axis.text.x = element_blank())

plot_CCLE_4_Multi_ARI <- plot_results_methods_combined(Results = CCLE_4_Multi_Res_ARI, NumGroups = 4, Title = "ARI (# true groups = 4)", Lower = -0.05, Upper = 0.7)
plot_CCLE_4_Multi_NMI <- plot_results_methods_combined(Results = CCLE_4_Multi_Res_NMI, NumGroups = 4, Title = "NMI (# true groups = 4)", Lower = -0.05, Upper = 0.7)
plot_CCLE_4_Multi_Geo <- plot_results_methods_combined(Results = CCLE_4_Multi_Res_Geo, NumGroups = 4, Title = NULL, Lower = -0.05, Upper = 0.7)
plot_CCLE_4_Multi_Geo <- plot_CCLE_4_Multi_Geo + theme(axis.text.x = element_blank())

ARI_NMI_4 <- ggarrange(plot_CCLE_4_Single_ARI, plot_CCLE_4_Single_NMI, plot_CCLE_4_Multi_ARI, plot_CCLE_4_Multi_NMI, ncol = 4, nrow = 1)
Geo_4 <- ggarrange(plot_CCLE_4_Single_Geo, plot_CCLE_4_Multi_Geo, ncol = 2, nrow = 1)

###########
#EINS_CCLE5
###########
#True labels
CCLE_5_True <- as.factor(EINS_CCLE5$Omics$Metadata$Proteomics$Site_Primary)
levels(CCLE_5_True) <- c(1:5)
CCLE_5_True <- as.numeric(CCLE_5_True)

#True order
CCLE_5_Sam <- rownames(EINS_CCLE5$Omics$Metadata$Proteomics)

#Single-omics data prep
CCLE_5_Single <- list()
CCLE_5_Ens <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE5$Ensemble$Samples$CHC$ClusterRes,
                               ClusterList = names(EINS_CCLE5$Ensemble$Samples$CHC$ClusterRes),
                               Labels = CCLE_5_True,
                               Samples = CCLE_5_Sam)
CCLE_5_Ens <- standardize_metrics(CCLE_5_Ens)
CCLE_5_Single$Ensemble <- CCLE_5_Ens

CCLE_5_Omics<- names(EINS_CCLE5$Single_Omics$HClustTree)
for(Omics in CCLE_5_Omics){
  Res <- ARI_NMI_HClust(Tree = EINS_CCLE5$Single_Omics$HClustTree[[Omics]],
                        Labels = CCLE_5_True,
                        Samples = CCLE_5_Sam)
  Res <- standardize_metrics(Res)
  CCLE_5_Single[[Omics]] <- Res
}

#Multi-omics data prep
CCLE_5_Multi <- list()
CCLE_5_Multi$Ensemble <- CCLE_5_Ens

CCLE_5_Methods <- names(EINS_CCLE5$Multi_Omics$ClusterRes)
for(Method in CCLE_5_Methods){
  if(Method == "SNF"){
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE5$Multi_Omics$ClusterRes$SNF$`euclidean squared`,
                            ClusterList = names(EINS_CCLE5$Multi_Omics$ClusterRes$SNF$`euclidean squared`),
                            Labels = CCLE_5_True,
                            Samples = CCLE_5_Sam)
    Res <- standardize_metrics(Res)
    CCLE_5_Multi[[Method]] <- Res
  }else if(Method == "GAUDI"){
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE5$Multi_Omics$ClusterRes$GAUDI$Kmeans,
                            ClusterList = names(EINS_CCLE5$Multi_Omics$ClusterRes$GAUDI$Kmeans),
                            Labels = CCLE_5_True,
                            Samples = CCLE_5_Sam)
    Res <- standardize_metrics(Res)
    CCLE_5_Multi[[Method]] <- Res
  }else {
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE5$Multi_Omics$ClusterRes[[Method]],
                            ClusterList = names(EINS_CCLE5$Multi_Omics$ClusterRes[[Method]]),
                            Labels = CCLE_5_True,
                            Samples = CCLE_5_Sam)
    Res <- standardize_metrics(Res)
    CCLE_5_Multi[[Method]] <- Res
  }
}
CCLE_5_Single <- CCLE_5_Single[c("Ensemble", "Methylation", "miRNA", "Transcriptomics", "Proteomics", "Metabolomics")]
CCLE_5_Multi <- CCLE_5_Multi[c("Ensemble", "MoCluster", "LRAcluster", "GAUDI", "COCA", "SNF", "MCIA", "MOFA")]

#combine results in single dataframe
CCLE_5_Single_Res <- results_combined(
  results_list = CCLE_5_Single,
  metrics = c("ARI", "NMI", "GeoMean"),
  line_sizes = c(2.8, 1, 1, 1, 1, 1)
)

CCLE_5_Single_Res_ARI <- CCLE_5_Single_Res %>% filter(Metric == "ARI")
CCLE_5_Single_Res_NMI <- CCLE_5_Single_Res %>% filter(Metric == "NMI")
CCLE_5_Single_Res_Geo <- CCLE_5_Single_Res %>% filter(Metric == "GeoMean")

CCLE_5_Multi_Res <- results_combined(
  results_list = CCLE_5_Multi,
  metrics = c("ARI", "NMI", "GeoMean"),
  facet_scales = "free_y",
  ncol = 2,
  palette = "Dark2",
  line_sizes = c(2.8, 1, 1, 1, 1, 1,1,1)
)


CCLE_5_Multi_Res_ARI <- CCLE_5_Multi_Res %>% filter(Metric == "ARI")
CCLE_5_Multi_Res_NMI <- CCLE_5_Multi_Res %>% filter(Metric == "NMI")
CCLE_5_Multi_Res_Geo <- CCLE_5_Multi_Res %>% filter(Metric == "GeoMean")

#plot results
plot_CCLE_5_Single_ARI <- plot_results_combined(Results = CCLE_5_Single_Res_ARI, NumGroups = 5, Title = "ARI (# true groups = 5)", Lower = -0.05, Upper = 0.6)
plot_CCLE_5_Single_NMI <- plot_results_combined(Results = CCLE_5_Single_Res_NMI, NumGroups = 5, Title = "NMI (# true groups = 5)", Lower = -0.05, Upper = 0.6)
plot_CCLE_5_Single_Geo <- plot_results_combined(Results = CCLE_5_Single_Res_Geo, NumGroups = 5, Title = NULL, Lower = -0.05, Upper = 0.6)
plot_CCLE_5_Single_Geo <- plot_CCLE_5_Single_Geo + theme(axis.text.x = element_blank())

plot_CCLE_5_Multi_ARI <- plot_results_methods_combined(Results = CCLE_5_Multi_Res_ARI, NumGroups = 5, Title = "ARI (# true groups = 5)", Lower = -0.05, Upper = 0.6)
plot_CCLE_5_Multi_NMI <- plot_results_methods_combined(Results = CCLE_5_Multi_Res_NMI, NumGroups = 5, Title = "NMI (# true groups = 5)", Lower = -0.05, Upper = 0.6)
plot_CCLE_5_Multi_Geo <- plot_results_methods_combined(Results = CCLE_5_Multi_Res_Geo, NumGroups = 5, Title = NULL, Lower = -0.05, Upper = 0.6)
plot_CCLE_5_Multi_Geo <- plot_CCLE_5_Multi_Geo + theme(axis.text.x = element_blank())

ARI_NMI_5 <- ggarrange(plot_CCLE_5_Single_ARI, plot_CCLE_5_Single_NMI, plot_CCLE_5_Multi_ARI, plot_CCLE_5_Multi_NMI, ncol = 4, nrow = 1)
Geo_5 <- ggarrange(plot_CCLE_5_Single_Geo, plot_CCLE_5_Multi_Geo, ncol = 2, nrow = 1)

###########
#EINS_CCLE6
###########
#True labels
CCLE_6_True <- as.factor(EINS_CCLE6$Omics$Metadata$Proteomics$Site_Primary)
levels(CCLE_6_True) <- c(1:6)
CCLE_6_True <- as.numeric(CCLE_6_True)

#True order
CCLE_6_Sam <- rownames(EINS_CCLE6$Omics$Metadata$Proteomics)

#Single-omics data prep
CCLE_6_Single <- list()
CCLE_6_Ens <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE6$Ensemble$Samples$CHC$ClusterRes,
                               ClusterList = names(EINS_CCLE6$Ensemble$Samples$CHC$ClusterRes),
                               Labels = CCLE_6_True,
                               Samples = CCLE_6_Sam)
CCLE_6_Ens <- standardize_metrics(CCLE_6_Ens)
CCLE_6_Single$Ensemble <- CCLE_6_Ens

CCLE_6_Omics<- names(EINS_CCLE6$Single_Omics$HClustTree)
for(Omics in CCLE_6_Omics){
  Res <- ARI_NMI_HClust(Tree = EINS_CCLE6$Single_Omics$HClustTree[[Omics]],
                        Labels = CCLE_6_True,
                        Samples = CCLE_6_Sam)
  Res <- standardize_metrics(Res)
  CCLE_6_Single[[Omics]] <- Res
}

#Multi-omics data prep
CCLE_6_Multi <- list()
CCLE_6_Multi$Ensemble <- CCLE_6_Ens

CCLE_6_Methods <- names(EINS_CCLE6$Multi_Omics$ClusterRes)
for(Method in CCLE_6_Methods){
  if(Method == "SNF"){
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE6$Multi_Omics$ClusterRes$SNF$`euclidean squared`,
                            ClusterList = names(EINS_CCLE6$Multi_Omics$ClusterRes$SNF$`euclidean squared`),
                            Labels = CCLE_6_True,
                            Samples = CCLE_6_Sam)
    Res <- standardize_metrics(Res)
    CCLE_6_Multi[[Method]] <- Res
  }else if(Method == "GAUDI"){
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE6$Multi_Omics$ClusterRes$GAUDI$Kmeans,
                            ClusterList = names(EINS_CCLE6$Multi_Omics$ClusterRes$GAUDI$Kmeans),
                            Labels = CCLE_6_True,
                            Samples = CCLE_6_Sam)
    Res <- standardize_metrics(Res)
    CCLE_6_Multi[[Method]] <- Res
  }else {
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE6$Multi_Omics$ClusterRes[[Method]],
                            ClusterList = names(EINS_CCLE6$Multi_Omics$ClusterRes[[Method]]),
                            Labels = CCLE_6_True,
                            Samples = CCLE_6_Sam)
    Res <- standardize_metrics(Res)
    CCLE_6_Multi[[Method]] <- Res
  }
}
CCLE_6_Single <- CCLE_6_Single[c("Ensemble", "Methylation", "miRNA", "Transcriptomics", "Proteomics", "Metabolomics")]
CCLE_6_Multi <- CCLE_6_Multi[c("Ensemble", "MoCluster", "LRAcluster", "GAUDI", "COCA", "SNF", "MCIA", "MOFA")]

#combine results in single dataframe
CCLE_6_Single_Res <- results_combined(
  results_list = CCLE_6_Single,
  metrics = c("ARI", "NMI", "GeoMean"),
  line_sizes = c(2.8, 1, 1, 1, 1, 1)
)

CCLE_6_Single_Res_ARI <- CCLE_6_Single_Res %>% filter(Metric == "ARI")
CCLE_6_Single_Res_NMI <- CCLE_6_Single_Res %>% filter(Metric == "NMI")
CCLE_6_Single_Res_Geo <- CCLE_6_Single_Res %>% filter(Metric == "GeoMean")

CCLE_6_Multi_Res <- results_combined(
  results_list = CCLE_6_Multi,
  metrics = c("ARI", "NMI", "GeoMean"),
  facet_scales = "free_y",
  ncol = 2,
  palette = "Dark2",
  line_sizes = c(2.8, 1, 1, 1, 1, 1,1,1)
)


CCLE_6_Multi_Res_ARI <- CCLE_6_Multi_Res %>% filter(Metric == "ARI")
CCLE_6_Multi_Res_NMI <- CCLE_6_Multi_Res %>% filter(Metric == "NMI")
CCLE_6_Multi_Res_Geo <- CCLE_6_Multi_Res %>% filter(Metric == "GeoMean")

#plot results
plot_CCLE_6_Single_ARI <- plot_results_combined(Results = CCLE_6_Single_Res_ARI, NumGroups = 6, Title = "ARI (# true groups = 6)", Lower = -0.05, Upper = 0.6)
plot_CCLE_6_Single_NMI <- plot_results_combined(Results = CCLE_6_Single_Res_NMI, NumGroups = 6, Title = "NMI (# true groups = 6)", Lower = -0.05, Upper = 0.6)
plot_CCLE_6_Single_Geo <- plot_results_combined(Results = CCLE_6_Single_Res_Geo, NumGroups = 6, Title = NULL, Lower = -0.05, Upper = 0.5)
plot_CCLE_6_Single_Geo <- plot_CCLE_6_Single_Geo + theme(axis.text.x = element_blank())

plot_CCLE_6_Multi_ARI <- plot_results_methods_combined(Results = CCLE_6_Multi_Res_ARI, NumGroups = 6, Title = "ARI (# true groups = 6)", Lower = -0.05, Upper = 0.6)
plot_CCLE_6_Multi_NMI <- plot_results_methods_combined(Results = CCLE_6_Multi_Res_NMI, NumGroups = 6, Title = "NMI (# true groups = 6)", Lower = -0.05, Upper = 0.6)
plot_CCLE_6_Multi_Geo <- plot_results_methods_combined(Results = CCLE_6_Multi_Res_Geo, NumGroups = 6, Title = NULL, Lower = -0.05, Upper = 0.5)
plot_CCLE_6_Multi_Geo <- plot_CCLE_6_Multi_Geo + theme(axis.text.x = element_blank())

ARI_NMI_6 <- ggarrange(plot_CCLE_6_Single_ARI, plot_CCLE_6_Single_NMI, plot_CCLE_6_Multi_ARI, plot_CCLE_6_Multi_NMI, ncol = 4, nrow = 1)
Geo_6 <- ggarrange(plot_CCLE_6_Single_Geo, plot_CCLE_6_Multi_Geo, ncol = 2, nrow = 1)

###########
#EINS_CCLE7
###########
#True labels
CCLE_7_True <- as.factor(EINS_CCLE7$Omics$Metadata$Proteomics$Site_Primary)
levels(CCLE_7_True) <- c(1:7)
CCLE_7_True <- as.numeric(CCLE_7_True)

#True order
CCLE_7_Sam <- rownames(EINS_CCLE7$Omics$Metadata$Proteomics)

#Single-omics data prep
CCLE_7_Single <- list()
CCLE_7_Ens <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE7$Ensemble$Samples$CHC$ClusterRes,
                               ClusterList = names(EINS_CCLE7$Ensemble$Samples$CHC$ClusterRes),
                               Labels = CCLE_7_True,
                               Samples = CCLE_7_Sam)
CCLE_7_Ens <- standardize_metrics(CCLE_7_Ens)
CCLE_7_Single$Ensemble <- CCLE_7_Ens

CCLE_7_Omics<- names(EINS_CCLE7$Single_Omics$HClustTree)
for(Omics in CCLE_7_Omics){
  Res <- ARI_NMI_HClust(Tree = EINS_CCLE7$Single_Omics$HClustTree[[Omics]],
                        Labels = CCLE_7_True,
                        Samples = CCLE_7_Sam)
  Res <- standardize_metrics(Res)
  CCLE_7_Single[[Omics]] <- Res
}

#Multi-omics data prep
CCLE_7_Multi <- list()
CCLE_7_Multi$Ensemble <- CCLE_7_Ens

CCLE_7_Methods <- names(EINS_CCLE7$Multi_Omics$ClusterRes)
for(Method in CCLE_7_Methods){
  if(Method == "SNF"){
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE7$Multi_Omics$ClusterRes$SNF$`euclidean squared`,
                            ClusterList = names(EINS_CCLE7$Multi_Omics$ClusterRes$SNF$`euclidean squared`),
                            Labels = CCLE_7_True,
                            Samples = CCLE_7_Sam)
    Res <- standardize_metrics(Res)
    CCLE_7_Multi[[Method]] <- Res
  }else if(Method == "GAUDI"){
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE7$Multi_Omics$ClusterRes$GAUDI$Kmeans,
                            ClusterList = names(EINS_CCLE7$Multi_Omics$ClusterRes$GAUDI$Kmeans),
                            Labels = CCLE_7_True,
                            Samples = CCLE_7_Sam)
    Res <- standardize_metrics(Res)
    CCLE_7_Multi[[Method]] <- Res
  }else {
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE7$Multi_Omics$ClusterRes[[Method]],
                            ClusterList = names(EINS_CCLE7$Multi_Omics$ClusterRes[[Method]]),
                            Labels = CCLE_7_True,
                            Samples = CCLE_7_Sam)
    Res <- standardize_metrics(Res)
    CCLE_7_Multi[[Method]] <- Res
  }
}
CCLE_7_Single <- CCLE_7_Single[c("Ensemble", "Methylation", "miRNA", "Transcriptomics", "Proteomics", "Metabolomics")]
CCLE_7_Multi <- CCLE_7_Multi[c("Ensemble", "MoCluster", "LRAcluster", "GAUDI", "COCA", "SNF", "MCIA", "MOFA")]

#combine results in single dataframe
CCLE_7_Single_Res <- results_combined(
  results_list = CCLE_7_Single,
  metrics = c("ARI", "NMI", "GeoMean"),
  line_sizes = c(2.8, 1, 1, 1, 1, 1)
)

CCLE_7_Single_Res_ARI <- CCLE_7_Single_Res %>% filter(Metric == "ARI")
CCLE_7_Single_Res_NMI <- CCLE_7_Single_Res %>% filter(Metric == "NMI")
CCLE_7_Single_Res_Geo <- CCLE_7_Single_Res %>% filter(Metric == "GeoMean")

CCLE_7_Multi_Res <- results_combined(
  results_list = CCLE_7_Multi,
  metrics = c("ARI", "NMI", "GeoMean"),
  facet_scales = "free_y",
  ncol = 2,
  palette = "Dark2",
  line_sizes = c(2.8, 1, 1, 1, 1, 1,1,1)
)


CCLE_7_Multi_Res_ARI <- CCLE_7_Multi_Res %>% filter(Metric == "ARI")
CCLE_7_Multi_Res_NMI <- CCLE_7_Multi_Res %>% filter(Metric == "NMI")
CCLE_7_Multi_Res_Geo <- CCLE_7_Multi_Res %>% filter(Metric == "GeoMean")

#plot results
plot_CCLE_7_Single_ARI <- plot_results_combined(Results = CCLE_7_Single_Res_ARI, NumGroups = 7, Title = "ARI (# true groups = 7)", Lower = -0.05, Upper = 0.55)
plot_CCLE_7_Single_NMI <- plot_results_combined(Results = CCLE_7_Single_Res_NMI, NumGroups = 7, Title = "NMI (# true groups = 7)", Lower = -0.05, Upper = 0.55)
plot_CCLE_7_Single_Geo <- plot_results_combined(Results = CCLE_7_Single_Res_Geo, NumGroups = 7, Title = NULL, Lower = -0.05, Upper = 0.5)
plot_CCLE_7_Single_Geo <- plot_CCLE_7_Single_Geo + theme(axis.text.x = element_blank())

plot_CCLE_7_Multi_ARI <- plot_results_methods_combined(Results = CCLE_7_Multi_Res_ARI, NumGroups = 7, Title = "ARI (# true groups = 7)", Lower = -0.05, Upper = 0.55)
plot_CCLE_7_Multi_NMI <- plot_results_methods_combined(Results = CCLE_7_Multi_Res_NMI, NumGroups = 7, Title = "NMI (# true groups = 7)", Lower = -0.05, Upper = 0.55)
plot_CCLE_7_Multi_Geo <- plot_results_methods_combined(Results = CCLE_7_Multi_Res_Geo, NumGroups = 7, Title = NULL, Lower = -0.05, Upper = 0.5)
plot_CCLE_7_Multi_Geo <- plot_CCLE_7_Multi_Geo + theme(axis.text.x = element_blank())

ARI_NMI_7 <- ggarrange(plot_CCLE_7_Single_ARI, plot_CCLE_7_Single_NMI, plot_CCLE_7_Multi_ARI, plot_CCLE_7_Multi_NMI, ncol = 4, nrow = 1)
Geo_7 <- ggarrange(plot_CCLE_7_Single_Geo, plot_CCLE_7_Multi_Geo, ncol = 2, nrow = 1)

###########
#EINS_CCLE8
###########
#True labels
CCLE_8_True <- as.factor(EINS_CCLE8$Omics$Metadata$Proteomics$Site_Primary)
levels(CCLE_8_True) <- c(1:8)
CCLE_8_True <- as.numeric(CCLE_8_True)

#True order
CCLE_8_Sam <- rownames(EINS_CCLE8$Omics$Metadata$Proteomics)

#Single-omics data prep
CCLE_8_Single <- list()
CCLE_8_Ens <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE8$Ensemble$Samples$CHC$ClusterRes,
                               ClusterList = names(EINS_CCLE8$Ensemble$Samples$CHC$ClusterRes),
                               Labels = CCLE_8_True,
                               Samples = CCLE_8_Sam)
CCLE_8_Ens <- standardize_metrics(CCLE_8_Ens)
CCLE_8_Single$Ensemble <- CCLE_8_Ens

CCLE_8_Omics<- names(EINS_CCLE8$Single_Omics$HClustTree)
for(Omics in CCLE_8_Omics){
  Res <- ARI_NMI_HClust(Tree = EINS_CCLE8$Single_Omics$HClustTree[[Omics]],
                        Labels = CCLE_8_True,
                        Samples = CCLE_8_Sam)
  Res <- standardize_metrics(Res)
  CCLE_8_Single[[Omics]] <- Res
}

#Multi-omics data prep
CCLE_8_Multi <- list()
CCLE_8_Multi$Ensemble <- CCLE_8_Ens

CCLE_8_Methods <- names(EINS_CCLE8$Multi_Omics$ClusterRes)
for(Method in CCLE_8_Methods){
  if(Method == "SNF"){
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE8$Multi_Omics$ClusterRes$SNF$`euclidean squared`,
                            ClusterList = names(EINS_CCLE8$Multi_Omics$ClusterRes$SNF$`euclidean squared`),
                            Labels = CCLE_8_True,
                            Samples = CCLE_8_Sam)
    Res <- standardize_metrics(Res)
    CCLE_8_Multi[[Method]] <- Res
  }else if(Method == "GAUDI"){
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE8$Multi_Omics$ClusterRes$GAUDI$Kmeans,
                            ClusterList = names(EINS_CCLE8$Multi_Omics$ClusterRes$GAUDI$Kmeans),
                            Labels = CCLE_8_True,
                            Samples = CCLE_8_Sam)
    Res <- standardize_metrics(Res)
    CCLE_8_Multi[[Method]] <- Res
  }else {
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE8$Multi_Omics$ClusterRes[[Method]],
                            ClusterList = names(EINS_CCLE8$Multi_Omics$ClusterRes[[Method]]),
                            Labels = CCLE_8_True,
                            Samples = CCLE_8_Sam)
    Res <- standardize_metrics(Res)
    CCLE_8_Multi[[Method]] <- Res
  }
}
CCLE_8_Single <- CCLE_8_Single[c("Ensemble", "Methylation", "miRNA", "Transcriptomics", "Proteomics", "Metabolomics")]
CCLE_8_Multi <- CCLE_8_Multi[c("Ensemble", "MoCluster", "LRAcluster", "GAUDI", "COCA", "SNF", "MCIA", "MOFA")]

#combine results in single dataframe
CCLE_8_Single_Res <- results_combined(
  results_list = CCLE_8_Single,
  metrics = c("ARI", "NMI", "GeoMean"),
  line_sizes = c(2.8, 1, 1, 1, 1, 1)
)

CCLE_8_Single_Res_ARI <- CCLE_8_Single_Res %>% filter(Metric == "ARI")
CCLE_8_Single_Res_NMI <- CCLE_8_Single_Res %>% filter(Metric == "NMI")
CCLE_8_Single_Res_Geo <- CCLE_8_Single_Res %>% filter(Metric == "GeoMean")

CCLE_8_Multi_Res <- results_combined(
  results_list = CCLE_8_Multi,
  metrics = c("ARI", "NMI", "GeoMean"),
  facet_scales = "free_y",
  ncol = 2,
  palette = "Dark2",
  line_sizes = c(2.8, 1, 1, 1, 1, 1,1,1)
)


CCLE_8_Multi_Res_ARI <- CCLE_8_Multi_Res %>% filter(Metric == "ARI")
CCLE_8_Multi_Res_NMI <- CCLE_8_Multi_Res %>% filter(Metric == "NMI")
CCLE_8_Multi_Res_Geo <- CCLE_8_Multi_Res %>% filter(Metric == "GeoMean")

#plot results
plot_CCLE_8_Single_ARI <- plot_results_combined(Results = CCLE_8_Single_Res_ARI, NumGroups = 8, Title = "ARI (# true groups = 8)", Lower = -0.05, Upper = 0.6)
plot_CCLE_8_Single_NMI <- plot_results_combined(Results = CCLE_8_Single_Res_NMI, NumGroups = 8, Title = "NMI (# true groups = 8)", Lower = -0.05, Upper = 0.6)
plot_CCLE_8_Single_Geo <- plot_results_combined(Results = CCLE_8_Single_Res_Geo, NumGroups = 8, Title = NULL, Lower = -0.05, Upper = 0.55)
plot_CCLE_8_Single_Geo <- plot_CCLE_8_Single_Geo + theme(axis.text.x = element_blank())

plot_CCLE_8_Multi_ARI <- plot_results_methods_combined(Results = CCLE_8_Multi_Res_ARI, NumGroups = 8, Title = "ARI (# true groups = 8)", Lower = -0.05, Upper = 0.6)
plot_CCLE_8_Multi_NMI <- plot_results_methods_combined(Results = CCLE_8_Multi_Res_NMI, NumGroups = 8, Title = "NMI (# true groups = 8)", Lower = -0.05, Upper = 0.6)
plot_CCLE_8_Multi_Geo <- plot_results_methods_combined(Results = CCLE_8_Multi_Res_Geo, NumGroups = 8, Title = NULL, Lower = -0.05, Upper = 0.55)
plot_CCLE_8_Multi_Geo <- plot_CCLE_8_Multi_Geo + theme(axis.text.x = element_blank())

ARI_NMI_8 <- ggarrange(plot_CCLE_8_Single_ARI, plot_CCLE_8_Single_NMI, plot_CCLE_8_Multi_ARI, plot_CCLE_8_Multi_NMI, ncol = 4, nrow = 1)
Geo_8 <- ggarrange(plot_CCLE_8_Single_Geo, plot_CCLE_8_Multi_Geo, ncol = 2, nrow = 1)

###########
#EINS_CCLE9
###########
#True labels
CCLE_9_True <- as.factor(EINS_CCLE9$Omics$Metadata$Proteomics$Site_Primary)
levels(CCLE_9_True) <- c(1:9)
CCLE_9_True <- as.numeric(CCLE_9_True)

#True order
CCLE_9_Sam <- rownames(EINS_CCLE9$Omics$Metadata$Proteomics)

#Single-omics data prep
CCLE_9_Single <- list()
CCLE_9_Ens <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE9$Ensemble$Samples$CHC$ClusterRes,
                               ClusterList = names(EINS_CCLE9$Ensemble$Samples$CHC$ClusterRes),
                               Labels = CCLE_9_True,
                               Samples = CCLE_9_Sam)
CCLE_9_Ens <- standardize_metrics(CCLE_9_Ens)
CCLE_9_Single$Ensemble <- CCLE_9_Ens

CCLE_9_Omics<- names(EINS_CCLE9$Single_Omics$HClustTree)
for(Omics in CCLE_9_Omics){
  Res <- ARI_NMI_HClust(Tree = EINS_CCLE9$Single_Omics$HClustTree[[Omics]],
                        Labels = CCLE_9_True,
                        Samples = CCLE_9_Sam)
  Res <- standardize_metrics(Res)
  CCLE_9_Single[[Omics]] <- Res
}

#Multi-omics data prep
CCLE_9_Multi <- list()
CCLE_9_Multi$Ensemble <- CCLE_9_Ens

CCLE_9_Methods <- names(EINS_CCLE9$Multi_Omics$ClusterRes)
for(Method in CCLE_9_Methods){
  if(Method == "SNF"){
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE9$Multi_Omics$ClusterRes$SNF$`euclidean squared`,
                            ClusterList = names(EINS_CCLE9$Multi_Omics$ClusterRes$SNF$`euclidean squared`),
                            Labels = CCLE_9_True,
                            Samples = CCLE_9_Sam)
    Res <- standardize_metrics(Res)
    CCLE_9_Multi[[Method]] <- Res
  }else if(Method == "GAUDI"){
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE9$Multi_Omics$ClusterRes$GAUDI$Kmeans,
                            ClusterList = names(EINS_CCLE9$Multi_Omics$ClusterRes$GAUDI$Kmeans),
                            Labels = CCLE_9_True,
                            Samples = CCLE_9_Sam)
    Res <- standardize_metrics(Res)
    CCLE_9_Multi[[Method]] <- Res
  }else {
    Res <- ARI_NMI_Clusters(ClusterRes = EINS_CCLE9$Multi_Omics$ClusterRes[[Method]],
                            ClusterList = names(EINS_CCLE9$Multi_Omics$ClusterRes[[Method]]),
                            Labels = CCLE_9_True,
                            Samples = CCLE_9_Sam)
    Res <- standardize_metrics(Res)
    CCLE_9_Multi[[Method]] <- Res
  }
}
CCLE_9_Single <- CCLE_9_Single[c("Ensemble", "Methylation", "miRNA", "Transcriptomics", "Proteomics", "Metabolomics")]
CCLE_9_Multi <- CCLE_9_Multi[c("Ensemble", "MoCluster", "LRAcluster", "GAUDI", "COCA", "SNF", "MCIA", "MOFA")]

#combine results in single dataframe
CCLE_9_Single_Res <- results_combined(
  results_list = CCLE_9_Single,
  metrics = c("ARI", "NMI", "GeoMean"),
  line_sizes = c(2.8, 1, 1, 1, 1, 1)
)

CCLE_9_Single_Res_ARI <- CCLE_9_Single_Res %>% filter(Metric == "ARI")
CCLE_9_Single_Res_NMI <- CCLE_9_Single_Res %>% filter(Metric == "NMI")
CCLE_9_Single_Res_Geo <- CCLE_9_Single_Res %>% filter(Metric == "GeoMean")

CCLE_9_Multi_Res <- results_combined(
  results_list = CCLE_9_Multi,
  metrics = c("ARI", "NMI", "GeoMean"),
  facet_scales = "free_y",
  ncol = 2,
  palette = "Dark2",
  line_sizes = c(2.8, 1, 1, 1, 1, 1,1,1)
)


CCLE_9_Multi_Res_ARI <- CCLE_9_Multi_Res %>% filter(Metric == "ARI")
CCLE_9_Multi_Res_NMI <- CCLE_9_Multi_Res %>% filter(Metric == "NMI")
CCLE_9_Multi_Res_Geo <- CCLE_9_Multi_Res %>% filter(Metric == "GeoMean")

#plot results
plot_CCLE_9_Single_ARI <- plot_results_combined(Results = CCLE_9_Single_Res_ARI, NumGroups = 9, Title = "ARI (# true groups = 9)", Lower = -0.05, Upper = 0.6)
plot_CCLE_9_Single_NMI <- plot_results_combined(Results = CCLE_9_Single_Res_NMI, NumGroups = 9, Title = "NMI (# true groups = 9)", Lower = -0.05, Upper = 0.6)
plot_CCLE_9_Single_Geo <- plot_results_combined(Results = CCLE_9_Single_Res_Geo, NumGroups = 9, Title = NULL, Lower = -0.05, Upper = 0.5)

plot_CCLE_9_Multi_ARI <- plot_results_methods_combined(Results = CCLE_9_Multi_Res_ARI, NumGroups = 9, Title = "ARI (# true groups = 9)", Lower = -0.05, Upper = 0.6)
plot_CCLE_9_Multi_NMI <- plot_results_methods_combined(Results = CCLE_9_Multi_Res_NMI, NumGroups = 9, Title = "NMI (# true groups = 9)", Lower = -0.05, Upper = 0.6)
plot_CCLE_9_Multi_Geo <- plot_results_methods_combined(Results = CCLE_9_Multi_Res_Geo, NumGroups = 9, Title = NULL, Lower = -0.05, Upper = 0.5)

ARI_NMI_9 <- ggarrange(plot_CCLE_9_Single_ARI, plot_CCLE_9_Single_NMI, plot_CCLE_9_Multi_ARI, plot_CCLE_9_Multi_NMI, ncol = 4, nrow = 1)
Geo_9 <- ggarrange(plot_CCLE_9_Single_Geo, plot_CCLE_9_Multi_Geo, ncol = 2, nrow = 1)


ARI_NMI <- ggarrange(ARI_NMI_3, ARI_NMI_4, ARI_NMI_5, ARI_NMI_6, ARI_NMI_7, ARI_NMI_8, ARI_NMI_9, nrow = 7)
ARI_NMI <- annotate_figure(ARI_NMI,
                bottom = "Number of clusters",
                left = "Metric value")


Geo <- ggarrange(Geo_3, Geo_4, Geo_5, Geo_6, Geo_7, Geo_8, Geo_9, nrow = 7)
Geo <- annotate_figure(Geo,
                       bottom = "Number of clusters",
                       left = "Geometric mean of ARI and NMI")