## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(EINS)
set.seed(47)

## -----------------------------------------------------------------------------
sapply(EINS_NCI60_Toy, dim)

## -----------------------------------------------------------------------------
EINS_NCI60 <- Ensemble_Integration_Suite$new()

## -----------------------------------------------------------------------------
EINS_NCI60$add_Omics_df(OmicsName = "Proteomics", DataDF = EINS_NCI60_Toy$Proteomics, MetaDF = EINS_NCI60_Toy$Metadata)
EINS_NCI60$add_Omics_df(OmicsName = "Transcriptomics", DataDF = EINS_NCI60_Toy$Transcriptomics, MetaDF = EINS_NCI60_Toy$Metadata)
EINS_NCI60$add_Omics_df(OmicsName = "Methylation", DataDF = EINS_NCI60_Toy$Methylation, MetaDF = EINS_NCI60_Toy$Metadata)

## -----------------------------------------------------------------------------
head(EINS_NCI60$Omics$Raw_Data$Proteomics)
head(EINS_NCI60$Omics$Metadata$Proteomics)

## -----------------------------------------------------------------------------
EINS_NCI60$pivot_Samples_Wide(SamplesInCol = TRUE)

## -----------------------------------------------------------------------------
EINS_NCI60$run_Preprocessing(OmicsName = "Proteomics", FunctionOrder = c("Normalization", "NAImpute"), NAMethod = "knn", NA_K_Neighbors = 10)
EINS_NCI60$run_Preprocessing(OmicsName = "Transcriptomics", FunctionOrder = c("NAImpute"), NAMethod = "knn", NA_K_Neighbors = 10)
EINS_NCI60$run_Preprocessing(OmicsName = "Methylation", FunctionOrder = c("Normalization", "NAImpute"), NAMethod = "knn", NA_K_Neighbors = 10)

## -----------------------------------------------------------------------------
head(EINS_NCI60$Preprocessed_Omics$Data$Proteomics)

## -----------------------------------------------------------------------------
EINS_NCI60$run_Sample_Matching()

## -----------------------------------------------------------------------------
head(EINS_NCI60$Preprocessed_Omics$Matched_Data$Proteomics)

## -----------------------------------------------------------------------------
EINS_NCI60$run_Single_Omics_Clustering(MinClust = 2, MaxClust = 8, HDModel = "all")

## -----------------------------------------------------------------------------
head(EINS_NCI60$Single_Omics$ClusterRes$Proteomics)

## -----------------------------------------------------------------------------
EINS_NCI60$run_Single_Omics_Hierarchical_Clustering()

## -----------------------------------------------------------------------------
head(EINS_NCI60$Single_Omics$DistanceMatrix$Proteomics)
plot(EINS_NCI60$Single_Omics$HClustTree$Proteomics)

## -----------------------------------------------------------------------------
EINS_NCI60$run_MoCluster(Clusters = 2:8, Components = 5, Linkage = "ward.D2")

## -----------------------------------------------------------------------------
head(EINS_NCI60$Multi_Omics$ClusterRes$MoCluster$Clusters_5)
head(EINS_NCI60$Multi_Omics$FeatureRes$MoCluster$Factors_5$Proteomics)

## -----------------------------------------------------------------------------
EINS_NCI60$run_MCIA(Components = 5, PerformCluster = TRUE, Clusters = 2:8)

## -----------------------------------------------------------------------------
head(EINS_NCI60$Multi_Omics$ClusterRes$MCIA$Clusters_5)
head(EINS_NCI60$Multi_Omics$FeatureRes$MCIA$Factors_5$Proteomics)

## -----------------------------------------------------------------------------
EINS_NCI60$run_LRAcluster(Clusters = 2:8, Linkage = "ward.D2")

## -----------------------------------------------------------------------------
head(EINS_NCI60$Multi_Omics$ClusterRes$LRAcluster$Clusters_5)

## -----------------------------------------------------------------------------
EINS_NCI60$run_Multi_Omics_Feature_HClust(nFactors = 5)

## -----------------------------------------------------------------------------
head(EINS_NCI60$Multi_Omics$Feature_HClust$DistanceMatrix$MoCluster)
plot(EINS_NCI60$Multi_Omics$Feature_HClust$HClustTree$MoCluster)

## -----------------------------------------------------------------------------
EINS_NCI60$run_Ensemble_Sample_MDS(Clusters = 5)

## -----------------------------------------------------------------------------
head(EINS_NCI60$Ensemble$Samples$MDS$ClusterRes$Clusters_5)
head(EINS_NCI60$Ensemble$Samples$MDS$DistanceMatrix$Clusters_5)
plot(EINS_NCI60$Ensemble$Samples$MDS$HClustTree$Clusters_5)

## -----------------------------------------------------------------------------
EINS_NCI60$run_Ensemble_Sample_CCA(Clusters = 5, Reference = "MoCluster", Factors = 5)

## -----------------------------------------------------------------------------
head(EINS_NCI60$Ensemble$Samples$CCA$Embedding_5$ClusterRes$Clusters_5)
head(EINS_NCI60$Ensemble$Samples$CCA$Embedding_5$EmbeddingMat)

## -----------------------------------------------------------------------------
EINS_NCI60$run_Ensemble_Integration_Feature(Method = "Average", nFactors = 5)

## -----------------------------------------------------------------------------
head(EINS_NCI60$Ensemble$Features$Average$DistanceMatrix)
plot(EINS_NCI60$Ensemble$Features$Average$HClustTree)

## -----------------------------------------------------------------------------
EINS_NCI60$run_Feature_Clustering(Clusters = 6, Method = "Average")

## -----------------------------------------------------------------------------
head(EINS_NCI60$Ensemble$Features$Average$ClusterRes)

## -----------------------------------------------------------------------------
EINS_NCI60$run_Over_Represenation(Method = "Average", OmicsName = "Proteomics", GeneNameType = "SYMBOL")

## -----------------------------------------------------------------------------
EINS_NCI60$Ensemble$Features$Average$ORA$Proteomics$Cluster_3$Enrichment

## -----------------------------------------------------------------------------
EINS_NCI60$plot_Heatmap_Cluster(FeatureClustering = TRUE, MethodResults = "MoCluster")

## -----------------------------------------------------------------------------
EINS_NCI60$Plots$Cluster_Heatmap$MoCluster$Clusters_5$Proteomics

## -----------------------------------------------------------------------------
EINS_NCI60$plot_Multi_Omics_Heatmap(MetadataFeatures = c("tissue of origin", "Epithelial", "sex", "p53"))

## -----------------------------------------------------------------------------
EINS_NCI60$Plots$Multi_Omics_Heatmap$Ensemble$MDS$Clusters_5

## -----------------------------------------------------------------------------
EINS_NCI60$plot_Scatterplot(XAxis = 1, YAxis = 2, MetadataFeature = "tissue of origin")

## -----------------------------------------------------------------------------
EINS_NCI60$Plots$Scatterplot$MoCluster$Clusters_5

## -----------------------------------------------------------------------------
EINS_NCI60$plot_Clusters_Metadata(MetadataFeatures = c("tissue of origin", "Epithelial", "sex", "p53"))

## -----------------------------------------------------------------------------
EINS_NCI60$Plots$Clusters_Metadata$Ensemble$MDS$Clusters_5

## -----------------------------------------------------------------------------
EINS_NCI60$plot_Sankey_Clusters(MetadataFeature = "tissue of origin")

## -----------------------------------------------------------------------------
EINS_NCI60$Plots$Sankey_Clusters$MoCluster

## -----------------------------------------------------------------------------
EINS_NCI60$plot_Sankey_Methods(MetadataFeature = "tissue of origin", Ensemble = TRUE, SNFDistance = "euclidean squared")

## -----------------------------------------------------------------------------
EINS_NCI60$Plots$Sankey_Methods$Selected$Clusters_5

## -----------------------------------------------------------------------------
EINS_NCI60$plot_Top_Feature_Weights(nFeatures = 10)

## -----------------------------------------------------------------------------
EINS_NCI60$Plots$Top_Feature_Weights$MoCluster$Factors_5$Proteomics

## -----------------------------------------------------------------------------
EINS_NCI60$plot_Enrich_Dot(Method = "Average", Categories = 5)

## -----------------------------------------------------------------------------
EINS_NCI60$Plots$ORA_Dot$Average$Proteomics$Cluster_3

## -----------------------------------------------------------------------------
EINS_NCI60$plot_Dendrogram(Data = "Single", MetadataHeight = c(-5, -5, -5), MetadataFeatures = c("tissue of origin", "Epithelial", "p53"), LabelMetadata = "tissue of origin", LabelSize = 0.5, BarLabelSize = 0.5)

## -----------------------------------------------------------------------------
EINS_NCI60$plot_Dendrogram(Data = "Ensemble MDS", MetadataHeight = -2, Clusters = 5, MetadataFeatures = c("tissue of origin", "Epithelial", "p53"), LabelMetadata = "tissue of origin", LabelSize = 0.5, BarLabelSize = 0.5)

## -----------------------------------------------------------------------------
EINS_NCI60$plot_Dendrogram(Data = "Ensemble Feature Average", MetadataHeight = -0.25)

## -----------------------------------------------------------------------------
EINS_NCI60$plot_Dendrogram(Data = "Multi-Omics Feature", MetadataHeight = -0.25)

