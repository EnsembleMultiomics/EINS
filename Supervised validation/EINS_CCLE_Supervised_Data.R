library(EINS)
set.seed(47)

#Set working directory with EINS CCLE data
setwd()
#Data import
CCLE_313_Meta <- read_csv("CCLE_313_Meta.csv")
CCLE_313_Prot <- read_csv("CCLE_313_Prot.csv")
CCLE_313_Metab <- read_csv("CCLE_313_Metab.csv")
CCLE_313_miRNA <- read_csv("CCLE_313_miRNA_Norm.csv")
CCLE_313_Meth <- read_csv("CCLE_313_Meth_Norm.csv")
CCLE_313_RNA <- read_csv("CCLE_313_RNA_Norm.csv")

CCLE_313_Prot <- column_to_rownames(CCLE_313_Prot, var = "Rows")
CCLE_313_Metab <- column_to_rownames(CCLE_313_Metab, var = "Rows")
CCLE_313_miRNA <- column_to_rownames(CCLE_313_miRNA, var = "Rows")
CCLE_313_Meth <- column_to_rownames(CCLE_313_Meth, var = "Rows")
CCLE_313_RNA <- column_to_rownames(CCLE_313_RNA, var = "Rows")
CCLE_313_Meta <- column_to_rownames(CCLE_313_Meta, var = "Rows")

##########
#For 9 groups only, for 17 groups use CCLE_313_ files
##########
#Group selection
##########
temp=(CCLE_313_Meta[,c(1,2,3,4)])
temp=temp[temp[,3]=="primary",]
temp=temp[temp[,4]%in%names(table(temp[,4]))[table(temp[,4])>=10],]
primary_site_temp=sort(unique(temp[,4]))

#Repeat for 3-9 groups, 3 groups given as example
CCLE_groups=rownames(temp[temp[,4]%in%primary_site_temp[1:9],])

CCLE_9_Meta <- CCLE_313_Meta[(rownames(CCLE_313_Meta) %in% CCLE_groups), ]
CCLE_9_Prot <- CCLE_313_Prot[,(colnames(CCLE_313_Prot) %in% CCLE_groups) ]
CCLE_9_Metab <- CCLE_313_Metab[,(colnames(CCLE_313_Metab) %in% CCLE_groups) ]
CCLE_9_miRNA <- CCLE_313_miRNA[,(colnames(CCLE_313_miRNA) %in% CCLE_groups) ]
CCLE_9_Meth <- CCLE_313_Meth[,(colnames(CCLE_313_Meth) %in% CCLE_groups) ]
CCLE_9_RNA <- CCLE_313_RNA[,(colnames(CCLE_313_RNA) %in% CCLE_groups) ]

##########
#END Group selection
##########

#Set wd for where to save RDS files
setwd()

#########
#EINS
#########
#Code only given for 17 groups, repeat in parallel session for 9 groups with CCLE_9_ files
#########
#Create R6 object for analysis 
EINS_CCLE <- Ensemble_Integration_Suite$new()

#Add raw data and metadata into R6 object
EINS_CCLE$add_Omics_df(OmicsName = "Proteomics", DataDF = CCLE_313_Prot, MetaDF = CCLE_313_Meta)
EINS_CCLE$add_Omics_df(OmicsName = "Metabolomics", DataDF = CCLE_313_Metab, MetaDF = CCLE_313_Meta)
EINS_CCLE$add_Omics_df(OmicsName = "miRNA", DataDF = CCLE_313_miRNA, MetaDF = CCLE_313_Meta)
EINS_CCLE$add_Omics_df(OmicsName = "Methylation", DataDF = CCLE_313_Meth, MetaDF = CCLE_313_Meta)
EINS_CCLE$add_Omics_df(OmicsName = "Transcriptomics", DataDF = CCLE_313_RNA, MetaDF = CCLE_313_Meta)

#----------------------Data Preparation----------------------
#Pivot matrices to ensure the samples are the columns (SamplesInCol = TRUE ==> in original data, samples are columns)
EINS_CCLE$pivot_Samples_Wide(SamplesInCol = TRUE)

#Preprocessing of the datasets
EINS_CCLE$run_Preprocessing(OmicsName = "Proteomics", FunctionOrder = c("Normalization", "NAImpute", "FilterCorrelation"), CorrelationThreshold = 0.9, NAMethod = "knn", NA_K_Neighbors = 10)
EINS_CCLE$run_Preprocessing(OmicsName = "Metabolomics", FunctionOrder = c("Normalization", "NAImpute", "FilterCorrelation"), CorrelationThreshold = 0.9, NAMethod = "knn", NA_K_Neighbors = 10)
EINS_CCLE$run_Preprocessing(OmicsName = "miRNA", FunctionOrder = c("Normalization", "NAImpute", "FilterCorrelation"), CorrelationThreshold = 0.9, NAMethod = "knn", NA_K_Neighbors = 10)
EINS_CCLE$run_Preprocessing(OmicsName = "Methylation", FunctionOrder = c("Normalization", "FilterCoverage", "NAImpute", "FilterCorrelation"), CorrelationThreshold = 0.9, CoverageThreshold = 0.8, NAMethod = "knn", NA_K_Neighbors = 10)
EINS_CCLE$run_Preprocessing(OmicsName = "Transcriptomics", FunctionOrder = c("Normalization", "FilterCoverage", "NAImpute", "FilterCorrelation"), CorrelationThreshold = 0.9, CoverageThreshold = 0.8, NAMethod = "knn", NA_K_Neighbors = 10)

#Matching samples across omics datasets
EINS_CCLE$run_Sample_Matching()

#----------------------Single omics clustering----------------------
#Single omics clustering with hierarchical clustering
EINS_CCLE$run_Single_Omics_Hierarchical_Clustering()

#----------------------Multi-omics integration----------------------
#Run MoCluster multi-omics integration
EINS_CCLE$run_MoCluster(Clusters = 2:20, Components = 9, Linkage = "ward.D2")

#Run MCIA multi-omics integration
EINS_CCLE$run_MCIA(PerformCluster = T, Clusters = 2:20, Components = 9)

#Run jNMF multi-omics integration
EINS_CCLE$run_jNMF(Clusters = 2:20, MaxIter = 200)

#Run iNMF multi-omics integration
EINS_CCLE$run_iNMF(Clusters = 2:20, MaxIter = 200)

#Run MOFA multi-omics integration
EINS_CCLE$run_MOFA(PerformCluster = T, Clusters = 2:20, Factors = 9)

#Run LRAcluster multi-omics integration
EINS_CCLE$run_LRAcluster(Clusters = 2:20, Linkage = "ward.D2")

#Run COCA multi-omics integration
EINS_CCLE$run_COCA(DatasetClusters = 9, FullClusters = 2:20, FullClusterMethod  = "hclust", FullHClustMethod = "ward.D2")

#Run SNF multi-omics integration
EINS_CCLE$run_SNF(Distance = c("euclidean squared"), Clusters = 2:20)

#Run GAUDI multi-omics integration with K-means clustering
EINS_CCLE$run_GAUDI(ClusterMethod = "Kmeans", KmeansClusters = 2:20, SingleOmicsUMAPDimensions = 15, ConcatUMAPDimensions = 9)

#----------------------Ensemble integration----------------------
#Run ensemble sample clustering for
for(i in 2:20){
  EINS_CCLE$run_Ensemble_Sample_CHC(Clusters = i)
}

saveRDS(EINS_CCLE, file = "CCLE_Super_Eval_17.rds")
