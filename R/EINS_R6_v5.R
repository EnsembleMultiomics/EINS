#' @title Ensemble_Integration_Suite Object.
#'
#' @description
#' An object in which to import multiple omics datasets, which can then be
#' analyzed with single omics, multi-omics integration and ensemble integration
#' methods. Visualization methods are included.
#'
#' @details
#' At least 2 different omics datasets are required, with metadata files. Single
#' omics analysis is possible with model-based clustering or hierarchical
#' clustering. 9 different multi-omics integration methods are available:
#' MoCluster, MCIA, jNMF, iNMF, LRAcluster, MOFA, COCA, GAUDI and SNF. After
#' multi-omics analysis, ensemble sample clustering and ensemble feature
#' clustering can be performed using the multi-omics integration results.
#' Visualization of the clustering results can then be performed.
#'
#' @export
Ensemble_Integration_Suite <- R6::R6Class(
  "Ensemble_Integration_Suite",
  public = list(
    #' @description
        #' Create a new Ensemble_Integration_Suite object.
    #' @return A new Ensemble_Integration_Suite object.
    initialize = function() {

    },

    #----------------------------------------------------------------ATTRIBUTES-------------------

    #Omics data list
    #' @field Omics List of different datasets per omics type.
    Omics = list(
      Metadata = NULL,
      Raw_Data = NULL,
      Sample_Data = NULL,
      QC_Data = NULL,
      Blank_Data = NULL,
      Pooled_Data = NULL
    ),

    #Adapted Omics data list
    #' @field Preprocessed_Omics List of datasets in proceeding steps of data
    #' preprocessing per omics type.
    Preprocessed_Omics = list(
      Data = NULL,
      Matched_Data = NULL,
      Non_Neg_Data = NULL
    ),

    #Outcome of data checks
    #' @field Data_Checks List of checks performed on omics data.
    Data_Checks = list(
      Optimal_Cluster_Number = NULL
    ),

    #Single-omics clustering results
    #' @field Single_Omics List of single omics analysis results per omics type.
    Single_Omics = list(
      Fit = NULL,
      ClusterRes = NULL,
      HClustTree = NULL,
      DistanceMatrix = NULL
    ),

    #Multi-omics integration list
    #' @field Multi_Omics List of multi-omics integration results per
    #' integration method.
    Multi_Omics = list(
      Fit = NULL,
      ClusterRes = NULL,
      FeatureRes = NULL,
      FactorRes = NULL,
      VarianceExp = NULL,
      CoordData = NULL
    ),

    #Ensemble integration list
    #' @field Ensemble List of ensemble integration results.
    Ensemble = list(
      Samples = NULL,
      Features = NULL,
      Factors = NULL
    ),

    #Plot list
    #' @field Plots List of created plots.
    Plots = list(
      Batch_PCA = NULL,
      Cluster_Heatmap = NULL,
      Multi_Omics_Heatmap = NULL,
      Scatterplot = NULL,
      Clusters_Metadata = NULL,
      Sankey_Clusters = NULL,
      Sankey_Methods = NULL,
      Top_Feature_Weights = NULL,
      Dendrogram = NULL,
      ORA_Dot = NULL
    ),

    #----------------------------------------------------------------DATA METHODS-----------------

    #Import methods

    ##Add omics per file
    #' @description
        #' Add individual omics raw data and metadata by  file path, to be
        #' repeated for each omics dataset.
    #' @param OmicsName Name of omics type.
    #' @param DataFile File path name for the raw omics data of the omics type.
    #' @param MetadataFile File path name for the metadata of the omics type.
    #' @param IndexColumnData Whether to use the first column as the index. Default
    #' is `TRUE`.
    #' @param IndexColumnMeta Whether to use the first column of the metadata
    #' file as the index. Default is `TRUE`.
    #' @returns Raw data stored in the object under $Omics$Raw_Data, metadata
    #' stored in the object under $Omics$Metadata
    add_Omics_file = function(OmicsName,
                              DataFile,
                              MetadataFile,
                              IndexColumnData = TRUE,
                              IndexColumnMeta = TRUE){
      meta_table = soda_read_table(file_path = MetadataFile, first_column_as_index = IndexColumnData)
      data_table = soda_read_table(file_path = DataFile, first_column_as_index = IndexColumnMeta)
      data_table = as.matrix(data_table)
      data_cols = colnames(data_table)
      meta_table = meta_table[match(data_cols, rownames(meta_table)),]

      self$Omics$Metadata[[OmicsName]] = meta_table
      self$Omics$Raw_Data[[OmicsName]] = data_table
    },

    ##Remove omics
    #' @description
        #' Remove all created datasets for a specific omics type. Removes both
        #' the raw datafiles associated with a specific omics type, as well as
        #' the preprocessed files. Analyses without this omics type will need to
        #' be run again.
    #' @param OmicsName Name of the omics type to be removed.
    #' @returns All $Omics and $Preprocessed_Omics files for the omics defined
    #' in `OmicsName` will become `NULL`, removing the matrices.
    remove_Omics = function(OmicsName){
      self$Omics$Metadata[[OmicsName]] = NULL
      self$Omics$Raw_Data[[OmicsName]] = NULL
      self$Omics$Sample_Data[[OmicsName]] = NULL
      self$Omics$QC_Data[[OmicsName]] = NULL
      self$Omics$Blank_Data[[OmicsName]] = NULL
      self$Omics$Pooled_Data[[OmicsName]] = NULL
      self$Preprocessed_Omics$Data[[OmicsName]] = NULL
      self$Preprocessed_Omics$Matched_Data[[OmicsName]] = NULL
      self$Preprocessed_Omics$Non_Neg_Data[[OmicsName]] = NULL
    },

    ##Add omics by df for building
    #' @description
        #' Add individual omics raw data and metadata from R environment, to be
        #' repeated for each omics type.
    #' @param OmicsName Name of omics type.
    #' @param DataDF Dataframe or matrix with raw omics data of the omics type.
    #' @param MetaDF Dataframe or matrix with metadata of the omics type.
    #' @returns Raw data stored in object under $Omics$Raw_Data, and metadata
    #' stored in object under $Omics$Metadata.
    add_Omics_df = function(OmicsName,
                            DataDF,
                            MetaDF){
      data_table = as.matrix(DataDF)
      meta_table = MetaDF
      data_cols = colnames(data_table)
      meta_table = meta_table[match(data_cols, rownames(meta_table)),]
      self$Omics$Raw_Data[[OmicsName]] = data_table
      self$Omics$Metadata[[OmicsName]] = meta_table
    },

    #Data check and modification methods
    #' @description
        #' Transpose the omics matrices. Required if input raw omics data does
        #' not store samples in columns and features in rows.
    #' @param SamplesInCol Are samples in columns? If FALSE (default), the raw
    #' omics matrices will be transposed.
    #' @returns $Omics$Raw_Data matrices with samples stored as columns.
    pivot_Samples_Wide = function(SamplesInCol){
      if(SamplesInCol == FALSE){ #samples are stored in rows
        self$Omics$Raw_Data = lapply(self$Omics$Raw_Data, t)
      }
    },


    #' @description
        #' Perform multiple different preprocessing and filtering steps on
        #' individual raw omics datasets. Implemented steps and order can be
        #' determined.
    #' @param OmicsName Name of omics type.
    #' @param FunctionOrder Character vector with any of the following 7 steps
    #' in any order.
    #'    * `"FeatureProp"`: Calculate the feature proportion of feature total
    #'    per sample.
    #'    * `"FilterCoverage"`: Filter features based on NA proportion across
    #'    samples, by user defined threshold.
    #'    * `"FilterBlankMean"`: Filter features based on its measurement in
    #'    blank samples. Multiply the mean feature measurement of the blank
    #'    samples by a user defined multiplier, and filter the features in which
    #'    less than the user defined threshold exceed the multiplied blank mean.
    #'    Only possible if blank samples are available.
    #'    * `"FilterCorrelation"`: Filter out correlated features. Calculate
    #'    pairwise correlations between features, and filter out the features
    #'    with correlation above the user defined threshold. One feature per
    #'    correlated group is kept in the dataset. The features which are
    #'    filtered out are stored in a list at
    #'    $Preprocessed_Omics$Correlated_Features.
    #'    * `"Normalization"`: Feature Z-score normalization.
    #'    * `"NAImpute"`: NA imputation or replacement with one of the following
    #'    methods: KNN imputation, missForest, zero replacement, omit NA.
    #'    * `"BatchCorrection"`: Correct for batch effect using ComBat. Only
    #'    possible if batch information is available for samples.
    #' @param QCName Common string in the name of quality control samples.
    #' @param BlankName Common string in the name of blank samples.
    #' @param PooledName Common string in the name of pooled samples.
    #' @param CoverageThreshold Numerical threshold for feature coverage
    #' filtering performed in `"FilterCoverage"`. Default is 0.8 (at least 0.8
    #' of the samples do not have NAs in a feature).
    #' @param BlankMeanMultiplier Numerical by which blank feature mean is
    #' multiplied in `"FilterBlankMean"`. Default is 2.
    #' @param BlankMeanThreshold Numerical threshold for blank mean filtering
    #' performed in `"FilterBlankMean"`. Default is 0.8 (at least 0.8 of the
    #' samples are larger than blank feature mean * BlankMeanMultiplier).
    #' @param CorrelationThreshold Numerical threshold for correlation filtering
    #' performed in `"FilterCorrelation"`. Default is 0.9 (if pairwise
    #' correlation between two features is 0.9 or larger, one of the features is
    #' removed).
    #' @param CorMean Whether to calculate the mean measurements for the
    #' correlated features. If `FALSE` (default), a representative feature from
    #' the correlated group of features is kept in the dataset. If `TRUE`, the
    #' mean of all feature measurements for a group of correlated features is
    #' included in the dataset.
    #' @param NAMethod String with method for `"NAImpute"`: `"KNN"`,
    #' `"missForest"`, `"zero replacement"` or `"omit NA"`.
    #' @param NA_K_Neighbors Numerical for the number of K nearest neighbors
    #' when using `"NAMethod"` `"KNN"`. Default is 10.
    #' @param BatchColumnName String of the column name which contains batching
    #' information in the metadata file.
    #' @param PreBatchPCA Whether to create a PCA plot of the sample data prior
    #' to batch effect correction. If `TRUE` (default), a PCA plot of the
    #' samples before batch effect correction will be created and printed.
    #' @param BatchPCAXAxis Numerical to indicate which principle component will
    #' be used for X axis in the PCA plot. Default is 1.
    #' @param BatchPCAYAxis Numerical to indicate which principle component will
    #' be used for Y axis in the PCA plot. Default is 2.
    #' @returns Preprocessed data matrix for the defined omics type. Stored in
    #' $Preprocessed_Omics$Data
    #' @references * Hastie T, Tibshirani R, Narasimhan B, Chu G (2025). impute:
    #'   impute: Imputation for microarray data. doi:10.18129/B9.bioc.impute,
    #'   R package version 1.82.0
    #'   * Stekhoven DJ (2022). missForest: Nonparametric Missing Value Imputation
    #'   using Random Forest. R package version 1.5.
    #'   * Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Zhang Y, Storey
    #'   JD, Torres LC (2025). sva: Surrogate Variable Analysis.
    #'   doi:10.18129/B9.bioc.sva, R package version 3.56.0
    run_Preprocessing = function(OmicsName, #name of omics dataset to process
                                 FunctionOrder = c("FeatureProp", "FilterCoverage", "FilterBlankMean", "FilterCorrelation", "Normalization", "NAImpute", "BatchCorrection"), #order of different subfunctions (remove subfunction name if not necessary)
                                 QCName = NULL, #Common part in name of quality control samples
                                 BlankName = NULL, #Common part in name of blank samples
                                 PooledName = NULL, #Common part in name of pooled samples
                                 CoverageThreshold = 0.8, #threshold for NA coverage
                                 BlankMeanMultiplier = 2, #multiplier for blank mean threshold
                                 BlankMeanThreshold = 0.8, #threshold for features below blank mean
                                 CorrelationThreshold = 0.9, #threshold for feature pairwise correlation
                                 CorMean = FALSE, #whether to calculate mean for correlated group or take representative feature
                                 NAMethod = NULL, #method for NA imputation (KNN, missForest, zero replacement, omit NA)
                                 NA_K_Neighbors = 10, #K nearest neighbors for KNN NA imputation
                                 BatchColumnName = NULL, #Name of column with batch information in metadata file
                                 PreBatchPCA = TRUE, #if pre batch correction PCA plot should be created and printed
                                 BatchPCAXAxis = 1, #X axis for Batch effect PCA plot
                                 BatchPCAYAxis = 2 #Y axis for Batch effect PCA plot
    ){
      #select the correct omics dataset
      Data = self$Omics$Raw_Data[[OmicsName]]
      Metadata = self$Omics$Metadata[[OmicsName]]

      #Create separate sample dataframes
      if(!is.null(QCName)){
        QCData = Data[, grep(QCName, colnames(Data))]
        self$Omics$QC_Data[[OmicsName]] = QCData
        Data = Data[, -grep(QCName, colnames(Data))]
      }
      if(!is.null(BlankName)){
        BlankData = Data[, grep(BlankName, colnames(Data))]
        self$Omics$Blank_Data[[OmicsName]] = BlankData
        Data = Data[, -grep(BlankName, colnames(Data))]
      }
      if(!is.null(PooledName)){
        PooledData = Data[, grep(PooledName, colnames(Data))]
        self$Omics$Pooled_Data[[OmicsName]] = PooledData
        Data = Data[, -grep(PooledName, colnames(Data))]
      }

      #store unprocessed versions of different sample dataframes
      self$Omics$Sample_Data[[OmicsName]] = Data

      #run through functions in order as specified in FunctionOrder
      for(Func in FunctionOrder){
        if(Func == "FeatureProp"){
          Data = Feature_Proportion(OmicSampleData = Data)
        }else if(Func == "FilterCoverage"){
          Data = filter_Coverage(OmicSampleData = Data,
                                 Coverage = CoverageThreshold)
        }else if(Func == "FilterBlankMean"){
          Data = filter_Blank_Mean(OmicSampleData = Data,
                                   OmicBlankData = BlankData,
                                   BlankMeanMultiplier = BlankMeanMultiplier,
                                   SampleThreshold = BlankMeanThreshold
          )
        }else if(Func == "Normalization"){
          Data = normalize_ZScore(OmicSampleData = Data)
        }else if(Func == "NAImpute"){
          Data = impute_NA(OmicSampleData = Data,
                           Method = NAMethod,
                           K_Neighbors = NA_K_Neighbors)
        }else if(Func == "BatchCorrection"){
          if(PreBatchPCA == TRUE){
            #PCA plot
            PCA_Pre_Title = paste0("PCA plot ", OmicsName, ", pre batch effect correction")
            tData = t(Data)
            PCA_Pre_Res = stats::prcomp(tData)
            PCA_Pre_Plot = PCA_Plot_Batch(Data = PCA_Pre_Res,
                                          Metadata = Metadata,
                                          BatchColumnName = BatchColumnName,
                                          XAxis = BatchPCAXAxis,
                                          YAxis = BatchPCAYAxis,
                                          Width = 5,
                                          Height = 5,
                                          YFontSize = 15,
                                          XFontSize = 15,
                                          Title = PCA_Pre_Title)
            self$Plots$Batch_PCA$Pre_Correction[[OmicsName]] = PCA_Pre_Plot
          }
          Data = Batch_Effect_Correction(OmicSampleData = Data,
                                         Metadata = Metadata,
                                         BatchColumn = BatchColumnName)
          #PCA plot
          PCA_Post_Title = paste0("PCA plot ", OmicsName, ", post batch effect correction")
          tData = t(Data)
          PCA_Post_Res = stats::prcomp(tData)
          PCA_Post_Plot = PCA_Plot_Batch(Data = PCA_Post_Res,
                                         Metadata = Metadata,
                                         BatchColumnName = BatchColumnName,
                                         XAxis = BatchPCAXAxis,
                                         YAxis = BatchPCAYAxis,
                                         Width = 5,
                                         Height = 5,
                                         YFontSize = 15,
                                         XFontSize = 15,
                                         Title = PCA_Post_Title)
          self$Plots$Batch_PCA$Post_Correction[[OmicsName]] = PCA_Post_Plot
        }else if(Func == "FilterCorrelation"){
          CorData = CorFilter(OmicSampleData = Data, CutOff = CorrelationThreshold, CorMean = CorMean)
          Data = CorData$RedMat
          self$Preprocessed_Omics$Correlated_Features[[OmicsName]] = CorData$CorGroups
        }
      }

      #save preprocessed data
      self$Preprocessed_Omics$Data[[OmicsName]] = Data
    },

    #Match samples
    #' @description
        #' Select the samples which match across the different omics datasets.
        #' Most methods require the same samples in the same order in all
        #' included omics datasets. This method also matches samples in the
        #' metadata.
    #' @returns If the samples do not match across all omics datasets, only the
    #' matching samples are stored in $Preprocessed_Omics$Matched_Data. The
    #' matrix columns are ordered the same way across omics datasets and
    #' metadata.
    #' @examples
        #' E <- Ensemble_Integration_Suite$new()
        #' E$add_Omics_df(OmicsName = "Proteomics",
        #' DataDF = EINS_NCI60_Toy$Proteomics, MetaDF = EINS_NCI60_Toy$Metadata)
        #' E$add_Omics_df(OmicsName = "Methylation",
        #' DataDF = EINS_NCI60_Toy$Methylation,
        #' MetaDF = EINS_NCI60_Toy$Methylation)=
        #' E$run_Sample_Matching()
        #' E$Preprocessed_Omics$Matched_Data$Proteomics
        #'
    run_Sample_Matching = function(){
      #Check for matched samples and save matched samples in own list
      if(!is.null(self$Preprocessed_Omics$Data)){
        Data = self$Preprocessed_Omics$Data
      }else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      OmicsNames = names(Data)
      Matched = list()
      SampleCounts = c()

      for(name in OmicsNames){
        Matched[[name]] = colnames(Data[[name]])
        SampleCounts = c(SampleCounts, ncol(Data[[name]]))
      }

      Matched = Reduce(intersect, Matched)
      MatchedCounts = length(Matched)

      if(MatchedCounts != max(SampleCounts)){
        for(name in OmicsNames){
          Matcheddf = Data[[name]][,Matched]
          self$Preprocessed_Omics$Matched_Data[[name]] = MatchedDF[,match(Matched, colnames(MatchedDF))]
        }
        warning(paste0("Different samples provided: ", length(Matched), " shared samples stored in Adapted_Omics_List$Matched_Data"))
      }else{
        for(name in OmicsNames){
          MatchedDF = Data[[name]]
          self$Preprocessed_Omics$Matched_Data[[name]] = MatchedDF[,match(Matched, colnames(MatchedDF))]
        }
      }
    },

    #optimal cluster number assessment with NbClust
    #' @description
        #' Calculate the optimal cluster number for integrated data using
        #' NbClust. Performs data integration with seven methods that provide
        #' sample by factor/component matrices: MoCluster, MCIA, jNMF, iNMF,
        #' LRAcluster, MOFA and GAUDI. Integration methods are performed with
        #' default parameters. A range of factors/components can be
        #' provided, for multiple integration results to be used. NbClust is
        #' then performed on these matrices. Multiple NbClust indexes are
        #' calculated per matrix. All results from NbClust indexes for all
        #' matrices are combined and the mean, mode and median optimal cluster
        #' number are calculated. The user can then determine the preferred
        #' cluster number.
    #' @param Components Numerical or range of the components/factors to perform
    #' the integrations with.
    #' @param NbClustDistance String with the distance measure to be used for
    #' the NbClust dissimilarity matrix computation. Must be one of:
    #' `"euclidean"`, `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"` or
    #' `"minkowski"`. Default is `"euclidean"`.
    #' @param NbClustMinClust Minimal number of clusters considered by NbClust.
    #' Value between 1 and (number of objects - 1). Default is 2.
    #' @param NbClustMaxClust Maximum number of clusters considered by NbClust.
    #' Value between 2 and (number of object -2), but larger than
    #' NbClustMinClust. Default is 12.
    #' @param NbClustMethod The cluster analysis method used by NbClust. Must be
    #' one of: `"ward.D"`, `"ward.D2"`, `"single"`, `"complete"`, `"average"`,
    #' `"mcquitty"`, `"median"`, `"centroid"`, `"kmeans"`. Default is
    #' `"ward.D2"`.
    #' @param NbClustIndex Indexes to be calculated by NbClust. Must be one of:
    #' `"kl"`, `"ch"`, `"hartigan"`, `"ccc"`, `"scott"`, `"marriot"`, `"trcovw"`,
    #' `"tracew"`, `"friedman"`, `"rubin"`, `"cindex"`, `"db"`, `"silhouette"`,
    #' `"duda"`, `"pseudot2"`, `"beale"`, `"ratkowsky"`, `"ball"`,
    #' `"ptbiserial"`, `"gap"`, `"frey"`, `"mcclain"`, `"gamma"`, `"gplus"`,
    #' `"tau"`, `"dunn"`, `"hubert"`, `"sdindex"`, `"dindex"`, `"sdbw"`, `"all"`
    #' (all indices except GAP, Gamma, Gplus and Tau), `"alllong"` (all indices
    #' with Gap, Gamma, Gplus and Tau included). Default is `"all"`.
    #' @param NbClustAlphaBeale Numerical, significance value for Beale's index.
    #' Default is 0.1.
    #' @returns A list of the optimal cluster number calculated by each included
    #' index for each integration method per component/factor. The mean, median
    #' and mode are calculated per integration method, and for all methods
    #' combined.
    #' @references   * Charrad M, Ghazzali N, Boiteau V, Niknafs A (2014).
    #'   “NbClust: An R Package for Determining the Relevant Number of Clusters
    #'   in a Data Set.” Journal of Statistical Software, 61(6), 1–36.
    #'   * Meng C (2025). mogsa: Multiple omics data integrative clustering and
    #'   gene set analysis. doi:10.18129/B9.bioc.mogsa, R package version 1.42.0
    #'   * Meng C, Kuster B, Culhane A, Gholami AM (2013). “A multivariate
    #'   approach to the integration of multi-omics datasets.” BMC Bioinformatics.
    #'   * Tsuyuzaki, K., and Nikaido, I. (2024). nnTensor: Non-Negative Tensor
    #'   Decomposition. R package version 1.3.0
    #'   * Chalise, P., Raghavan, R., and Fridley, B. (2025). IntNMF: Integrative
    #'   Clustering of Multiple Genomic Dataset. R package version 1.3.0
    #'   * Lu, X., Meng, J., Zhou, Y., Jiang, L., and Yan, F. (2020).MOVICS: an R
    #'   package for multi-omics integration and visualization in cancer
    #'   subtyping. Bioinformatics, btaa1018.
    #'   * Wu D, Wang D, Zhang MQ, Gu J (2015). Fast dimension reduction and
    #'   integrative clustering of multi-omics data using low-rank approximation:
    #'   application to cancer molecular classification. BMC Genomics, 16(1):1022.
    #'   * Argelaguet R, Velten B, Arnol D, Dietrich S, Zenz T, Marioni JC,
    #'   Buettner F, Huber W, Stegle O (2018). “Multi‐Omics Factor Analysis—a
    #'   framework for unsupervised integration of multi‐omics data sets.”
    #'   Molecular Systems Biology, 14.
    #'   * Castellano-Escuder P, Zachman DK, Han K, Hirschey MD. GAUDI:
    #'   interpretable multi-omics integration with UMAP embeddings and
    #'   density-based clustering. Nat Commun. 2025 Jul 1;16(1):5771.
    run_Optimal_Cluster_Number = function(Components = NULL, #single value or range of possible components
                                          NbClustDistance = "euclidean",
                                          NbClustMinClust = 2,
                                          NbClustMaxClust = 12,
                                          NbClustMethod = "ward.D2",
                                          NbClustIndex = "all",
                                          NbClustAlphaBeale = 0.1){
      #Data checks
      if(!is.null(self$Preprocessed_Omics$Matched_Data)){
        Data = self$Preprocessed_Omics$Matched_Data
      } else if(!is.null(self$Preprocessed_Omics$Data)){
        Data = self$Preprocessed_Omics$Data
      } else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      OptimalClusterNumber = Optimal_Cluster_Number(Data = Data,
                                                    Components = Components,
                                                    NbClustDistance = NbClustDistance,
                                                    NbClustMinClust = NbClustMinClust,
                                                    NbClustMaxClust = NbClustMaxClust,
                                                    NbClustMethod = NbClustMethod,
                                                    NbClustIndex = NbClustIndex,
                                                    NbClustAlphaBeale = NbClustAlphaBeale)

      #save results
      self$Data_Checks$Optimal_Cluster_Number = OptimalClusterNumber
    },

    #' @description
        #' Uses model-based clustering method HDDC from the HDclassif package
        #' to provide clustering results for individual omics datasets.
        #' This is a specific parametrization of the Gaussian mixture model to
        #' work with high dimensional data. When running
        #' `run_Single_Omics_Clustering`, clustering is performed on all omics
        #' datasets imported into the object.
    #' @param MinClust Minimum number of clusters considered by the model.
    #' Default is 2.
    #' @param MaxClust Maximum number of clusters considered by the model.
    #' Value needs to be at least (MinClust + 1). Default is 12.
    #' @param HDModel Name of the model. Must be one or combination of:
    #' `"AkjBkQkDk"` (default), `"AkBkQkDk"`, `"ABkQkDk"`, `"AkjBQkDk"`,
    #' `"AkBQkDk"`, `"ABQkDk"`, `"AkjBkQkD"`, `"AkBkQkD"`, `"ABkQkD"`,
    #' `"AkjBQkD"`, `"AkBQkD"`, `"ABQkD"`, `"AjBQD"` or `"ABQD"`. If a
    #' combination of models is selected, the model which maximizes the BIC
    #' criterion is kept. All models can be run with `"ALL"`.
    #' @param HDThreshold Threshold used in Cattell's Scree-Test. Must be
    #' between 0 and 1. Default is 0.2.
    #' @param HDCriterion Criterion used to select the best model. Either
    #' `"BIC"` (default) or `"ICL"`.
    #' @param HDItermax Maximum number of iteration allowed. Default is 200.
    #' @param HDEps Double positive, default is 0.001. Used as the stopping
    #' criterion: the algorithm stops if the difference between two successive
    #' log-likelihoods is lower than `"HDEps"`.
    #' @param HDAlgo The algorithm to be used. Must be one of: `"EM"`
    #' (expectation-maximization), `"CEM"` (classification EM), or `"SEM"`
    #' (stochastic EM). Default is `"EM"`.
    #' @param HDD_select Method to select the intrinsic dimension. Must be one
    #' of: `"Cattell"` (default), or `"BIC"`.
    #' @param HDInit Method to initialize the EM algorithm. Must be one of:
    #' `"kmeans"` (default), `"param"`, `"random"`, `"mini-em"` or `"vector"`.
    #' If `"vector"` is selected, argument `"HDInit.vector"` should be added.
    #' @param HDInit.vector Vector of integers or factors of the same length as
    #' the number of samples. Only used if `"HDInit"` is `"vector"`.
    #' @param HDMini.nb Vector of 2 integers. Parameter used in `"mini-em"`
    #' initialization. First integer indicates the amount of times the algorithm
    #' will be launched, the second integer indicates how many iterations are
    #' performed each launch. Default is `"c(5,10)"`.
    #' @param HDScaling Whether to perform scaling on the dataset. Default is
    #' `FALSE`.
    #' @param HDMin.individuals The minimum number of samples per cluster. Must
    #' be 2 (default) or greater.
    #' @param HDNoise.ctrl Parameter to minimize the intrinsic dimension
    #' selected. Dimensions with eigenvalues lower than `"HDNoise.ctrl"` are not
    #' used. Default is 1e-8.
    #' @param HDMc.cores Whether parallel computing is used. Default is 1,
    #' indicating no parallel computing. If `"HDMc.cores"` is larger than 1,
    #' parallel computing is used with `"HDMc.cores"` cores.
    #' @param HDNb.rep Number of times each estimation is repeated, with the
    #' estimation with the highest log-likelihood kept. Default is 1 (no
    #' repetition).
    #' @param HDKeepAllRes Whether results of all runs should be kept. Default
    #' is `TRUE`.
    #' @param HDKmeans.control List of kmeans initialization parameters:
    #' `"iter.max"`, `"nstart"` and `"algorithm"`.
    #' @param HDD_max Maximum number of dimensions to be computed. Default is
    #' 100.
    #' @param HDSubset Clustering can be performed on a subset of the data. If
    #' `"HDSubset"` is smaller than the number of samples, clustering is
    #' performed on a random subset of that size. Then, the posterior of the
    #' clustering is computed on the full sample. Default is `Inf`, indicating
    #' no subsetting.
    #' @returns The arguments and fit of the model for each omics dataset,
    #' stored in $Single_Omics$Fit. A dataframe with the cluster assignment for
    #' each omics dataset, stored in $Single_Omics$ClusterRes.
    #' @references Bergé L, Bouveyron C, Girard S (2012). “HDclassif: An R
    #' Package for Model-Based Clustering and Discriminant Analysis of
    #' High-Dimensional Data.” Journal of Statistical Software, 46(6), 1–29.
    run_Single_Omics_Clustering = function(MinClust = 2,
                                           MaxClust = 12,
                                           HDModel = "AkjBkQkDk",
                                           HDThreshold = 0.2,
                                           HDCriterion = "BIC",
                                           HDItermax = 200,
                                           HDEps = 0.001,
                                           HDAlgo = "EM",
                                           HDD_select = "Cattell",
                                           HDInit = "kmeans",
                                           HDInit.vector,
                                           HDMini.nb = c(5, 10),
                                           HDScaling = FALSE,
                                           HDMin.individuals = 2,
                                           HDNoise.ctrl = 1e-08,
                                           HDMc.cores = 1,
                                           HDNb.rep = 1,
                                           HDKeepAllRes = TRUE,
                                           HDKmeans.control = list(),
                                           HDD_max = 100,
                                           HDSubset = Inf
    ){
      #data selection
      #No matching needed, so unmatched data used if available
      if(!is.null(self$Preprocessed_Omics$Data)){
        Data = self$Preprocessed_Omics$Data
      } else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      #HDclassif requires samples in row
      tData = lapply(Data, t)

      OmicsNames = names(tData)
      for(omics in OmicsNames){
        SingleData = tData[[omics]]
        SingleDataClustering = HDclassif::hddc(data = SingleData,
                                               K = MinClust:MaxClust,
                                               model = HDModel,
                                               threshold = HDThreshold,
                                               criterion = HDCriterion,
                                               itermax = HDItermax,
                                               eps = HDEps,
                                               algo = HDAlgo,
                                               d_select = HDD_select,
                                               init = HDInit,
                                               init.vector = HDInit.vector,
                                               mini.nb = HDMini.nb,
                                               scaling = HDScaling,
                                               min.individuals = HDMin.individuals,
                                               noise.ctrl = HDNoise.ctrl,
                                               mc.cores = HDMc.cores,
                                               nb.rep = HDNb.rep,
                                               keepAllRes = HDKeepAllRes,
                                               kmeans.control = HDKmeans.control,
                                               d_max = HDD_max,
                                               subset = HDSubset
        )
        self$Single_Omics$Fit[[omics]] = SingleDataClustering
        self$Single_Omics$ClusterRes[[omics]] = data.frame(row.names = rownames(SingleData),
                                                                Cluster = SingleDataClustering$class,
                                                                stringsAsFactors = FALSE)
      }
    },

    #' @description
        #' Sample hierarchical clustering to be performed on each individual
        #' omics dataset. When running `run_Single_Omics_Hierarchical_Clustering`,
        #' hierarchical clustering is performed on all omics datasets imported
        #' into the object.
    #' @param Distance Distance measure to be used. Must be one of:
    #' `"euclidean"` (default), `"maximum"`, `"manhattan"`, `"canberra"`,
    #' `"binary"` or `"minkowski"`. If `"minkowski"` is selected, argument
    #' `"MinkowskiPower"` needs to be included.
    #' @param MinkowskiPower Power of the Minkowski distance. Default is NULL.
    #' @param Linkage Agglomeration method to be used. Must be one of: `"ward.D"`,
    #' `"ward.D2"` (default), `"single"`, `"complete"`, `"average"`, `"mcquitty"`,
    #' `"median"` or `"centroid"`.
    #' @returns Distance matrix for each individual omics dataset, stored in
    #' $Single_Omics$DistanceMatrix. Hclust object for each individual omics
    #' dataset, stored in $Single_Omics$HClustTree. Dendrogram for each
    #' individual omics dataset, stored in $Plots$Dendrogram$Single_Omics.
    run_Single_Omics_Hierarchical_Clustering = function(Distance = "euclidean",
                                                        MinkowskiPower = NULL,
                                                        Linkage = "ward.D2"){
      #data selection
      #No matching needed, so unmatched data used if available
      if(!is.null(self$Preprocessed_Omics$Data)){
        Data = self$Preprocessed_Omics$Data
      } else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      #dist creates distance matrix for rows, so need samples in rows
      tData = lapply(Data, t)

      OmicsName = names(tData)
      for(omics in OmicsName){
        SingleData = tData[[omics]]
        SingleOmicsDist = stats::dist(SingleData,
                                      method = Distance,
                                      p = MinkowskiPower)
        SingleOmicsHClust = stats::hclust(SingleOmicsDist,
                                          method = Linkage)
        self$Single_Omics$HClustTree[[omics]] = SingleOmicsHClust
        self$Single_Omics$DistanceMatrix[[omics]] = SingleOmicsDist
        self$Plots$Dendrogram$Single_Omic[[omics]] = stats::as.dendrogram(SingleOmicsHClust)
      }
    },

    #----------------------------------------------------------------INTEGRATION METHODS----------

    #MoCluster
    #' @description
        #' Perform multi-omics integration and sample clustering using MoCluster
        #' from the mogsa package.
    #' @param Clusters The number of sample clusters to be created. Can be a
    #' single integer to calculate a single clustering assignment, or a range to
    #' calculate multiple clustering assignment.
    #' @param Components The number of components to reduce the data to. Needs
    #' be a single integer.
    #' @param Method String to indicate the integration approach. Must be one of:
    #' `"globalScore"` (consensus PCA), `"blockScore"` (generalized canonical
    #' correlation analysis) or `"blockLoading"` (multiple co-inertia analysis).
    #' Default is `"globalScore"`.
    #' @param Option String to indicate how matrices should be normalized. Must
    #' be one of: `"lambda1"` (matrix divided by first singular value),
    #' `"inertia"` (matrix divided by its total inertia) or `"uniform"` (no
    #' normalization). Default is `"lambda1"`.
    #' @param k The number of non-zero coefficients for the variable loading
    #' vectors to determine the sparsity. If `k` >= 1, `k` is the absolute
    #' number of non-zero coefficients. If 0  < `k` < 1, `k` is the proportion
    #' of non-zero coefficients. A single value or vector with length equal to
    #' the number of omics datasets can be provided. Default is `"all"`,
    #' indicating no sparsity.
    #' @param Center Whether the variables should be centered. Default is `TRUE`.
    #' @param Scale Whether the variables should be scaled. Default is `FALSE`.
    #' @param MaxIter Maximum number of iterations in the algorithm. Default is
    #' 1000.
    #' @param svdSolver String to indicate the method used for singular value
    #' decomposition. Must be one of: `"fast.svd"` (default), `"svd"` or
    #' `"propack"`. `"fast.svd"` has a good compromise between speed and
    #' robustness compared to `"svd"` (slower but more robust) and `"propack"`
    #' (faster but less likely to converge).
    #' @param Distance Distance measure to be used for distance matrix
    #' calculation of the integrated data. Must be one of: `"euclidean"`
    #' (default), `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"` or
    #' `"minkowski"`. If `"minkowski"` is selected, argument `"MinkowskiPower"`
    #' needs to be included.
    #' @param MinkowskiPower Power of the Minkowski distance. Default is NULL.
    #' @param Linkage Agglomeration method to be used for the hierarchical
    #' clustering. Must be one of: `"ward.D"`, `"ward.D2"`, `"single"`,
    #' `"complete"` (default), `"average"`, `"mcquitty"`, `"median"` or
    #' `"centroid"`.
    #' @returns A list of results. Fit of the integration for each included
    #' cluster number, stored in $Multi_Omics$Fit$MoCluster. Sample cluster
    #' assignment for each included cluster number, stored in
    #' $Multi_Omics$ClusterRes$MoCluster. Loading vectors (feature weights)
    #' calculated for each omics dataset, stored by component number in
    #' $Multi_Omics$FeatureRes$MoCluster. Variance in the sample clusters
    #' explained by each component for each included cluster number, stored in
    #' $Multi_Omics$VarianceExp$MoCluster. Sample factor score for each
    #' included cluster number, stored in $Multi_Omics$CoordData$MoCluster.
    #' @references Meng C (2025). mogsa: Multiple omics data integrative
    #' clustering and gene set analysis. doi:10.18129/B9.bioc.mogsa, R package
    #' version 1.42.0
    run_MoCluster = function(Clusters = NULL, #number of clusters
                             Components = NULL, #number of latent variables to calculate
                             Method = "globalScore", #calculation methods: globalScore = cPCA (default), blockScore = gCCA, blockLoading = MCIA
                             Option = "lambda1", #option for data normalization: lambda1 (default), inertia, uniform
                             k = "all", #absolute number or proportion of non-zero coefficients for loading vectors, default = all (no sparsity)
                             Center = TRUE, #logical to indicate centering
                             Scale = FALSE, #logical to indicate scaling
                             MaxIter = 1000,
                             svdSolver = "fast.svd", #fast.svd (default), svd (slowest) or propack (fastest)
                             Distance = "euclidean", #distance matrix distance
                             MinkowskiPower = NULL,
                             Linkage = "complete" #hclust agglomeration method: ward.D, ward.D2, single, complete (default), average, mcquitty, median, centroid
    ){
      #Data checks
      ##Same number of samples in columns needed of MoCluster and no NA allowed
      #first, check for matched and no NA, then matched only, then no NA only
      if(!is.null(self$Preprocessed_Omics$Matched_Data)){
        Data = self$Preprocessed_Omics$Matched_Data
      } else if(!is.null(self$Preprocessed_Omics$Data)){
        Data = self$Preprocessed_Omics$Data
      } else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      #run mogsa steps
      moas = mogsa::mbpca(x = Data,
                          ncomp = Components,
                          k = k,
                          method = Method,
                          option = Option,
                          center = Center,
                          scale = Scale,
                          moa = TRUE,
                          maxiter = MaxIter,
                          svd.solver = svdSolver,
                          verbose = FALSE)
      scrs = mogsa::moaScore(moa = moas)
      MCdist = stats::dist(x = scrs,
                           method = Distance,
                           p = MinkowskiPower)
      MCdend = stats::hclust(MCdist,
                             method = Linkage)

      #save feature results
      MCfeatres = moas@loading[which(rowSums(moas@loading != 0) > 0), ]
      colnames(MCfeatres)[1:ncol(MCfeatres)] <- paste0("V", 1:ncol(MCfeatres))
      MCf = sub('_[^_]*$', '', rownames(MCfeatres))
      MCd = sub('.*_', '', rownames(MCfeatres))
      rownames(MCfeatres) = NULL
      MCFeatFull = as.data.frame(cbind(Features = MCf, Dataset = MCd,MCfeatres))
      MCOmicsNames = names(Data)
      for(name in MCOmicsNames){
        MCFeatName = subset(MCFeatFull, Dataset == name)
        self$Multi_Omics$FeatureRes$MoCluster[[paste0("Factors_", Components)]][[name]] = subset(MCFeatName, select = -Dataset)
      }

      #save explained variance
      MCvarex <- as.data.frame(moas@ctr.tab)
      colnames(MCvarex)[1:ncol(MCvarex)] = paste0("V", 1:ncol(MCvarex))
      self$Multi_Omics$VarianceExp$MoCluster[[paste0("Factors_", Components)]] = MCvarex
      Coord = as.data.frame(moas@fac.scr)
      rownames(Coord) = colnames(Data[[1]])
      self$Multi_Omics$CoordData$MoCluster[[paste0("Factors_", Components)]] = Coord

      #save results
      for(i in Clusters){
        ClusterRes = data.frame(row.names = colnames(Data[[1]]),
                                Cluster = stats::cutree(MCdend,
                                                        k = i),
                                stringsAsFactors = FALSE)
        self$Multi_Omics$ClusterRes$MoCluster[[paste0("Clusters_", i)]] = ClusterRes
        self$Multi_Omics$Fit$MoCluster[[paste0("Clusters_", i)]] = moas
      }
    },

    #MCIA
    #' @description
        #' Perform multi-omics integration using multiple co-inertia analysis
        #' (MCIA) from the omicade4 package.
    #' @param PerformCluster Whether to perform k-means clustering on the sample
    #' embeddings. Default is `TRUE`.
    #' @param Clusters The number of sample clusters to be created. Can be a
    #' single integer to calculate a single cluster assignment, or a range to
    #' calculate multiple cluster assignments.
    #' @param Components The number of components to reduce the data to. Needs
    #' to be a single integer.
    #' @param SVD Whether singular value decomposition should be used. Default
    #' is `TRUE`.
    #' @returns A list of results. Fit of the integration, stored by component
    #' number in $Multi_Omics$Fit$MCIA. If PerformCluster = `TRUE`, sample cluster
    #' assignment for each included cluster number, stored in
    #' $Multi_Omics$ClusterRes$MCIA Loading vectors (feature weights)
    #' calculated for each omics dataset, stored by component number in
    #' $Multi_Omics$FeatureRes$MCIA. Sample factor score, stored by component
    #' number in $Multi_Omics$CoordData$MCIA.
    #' @references Meng C, Kuster B, Culhane A, Gholami AM (2013). “A
    #' multivariate approach to the integration of multi-omics datasets.” BMC
    #' Bioinformatics.
    run_MCIA = function(PerformCluster = TRUE, #whether to perform clustering
                        Clusters = NULL, #number of clusters
                        Components = NULL, #number of axes to keep
                        SVD = TRUE #whether to use svd
    ){

      #data checks
      ##Same number of samples in columns needed for MCIA
      if(!is.null(self$Preprocessed_Omics$Matched_Data)){
        Data = self$Preprocessed_Omics$Matched_Data
      } else if(!is.null(self$Preprocessed_Omics$Data)){
        Data = self$Preprocessed_Omics$Data
      } else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      #run MCIA
      mcoin <- omicade4::mcia(df.list = Data,
                              cia.nf = Components,
                              cia.scan = FALSE,
                              nsc = TRUE,
                              svd = SVD)

      #no clustering results
      #save feature results
      MCIAOmicsNames = names(Data)
      MCIAfeatres = mcoin$mcoa$Tco[which(rowSums(mcoin$mcoa$Tco != 0) > 0), ]
      colnames(MCIAfeatres)[1:ncol(MCIAfeatres)] = paste0("V", 1:ncol(MCIAfeatres))
      rownames(MCIAfeatres) = NULL
      MCIAd = mcoin$mcoa$TC$"T"
      FeatNames = c()
      for(omics in MCIAOmicsNames){
        FeatNames = c(FeatNames, rownames(Data[[omics]]))
      }
      MCIAf = FeatNames
      MCIAfeatres <- cbind(Features = MCIAf,
                           Dataset = MCIAd,
                           MCIAfeatres)

      for(name in MCIAOmicsNames){
        MCIAFeatName = subset(MCIAfeatres, Dataset == name)
        self$Multi_Omics$FeatureRes$MCIA[[paste0("Factors_", Components)]][[name]] = subset(MCIAFeatName, select = -Dataset)
      }
      Coord = as.data.frame(mcoin$mcoa$SynVar)
      rownames(Coord) = colnames(Data[[1]])
      self$Multi_Omics$CoordData$MCIA[[paste0("Factors_", Components)]] = Coord

      #create clusters and save results
      if(PerformCluster == TRUE){
        for(i in Clusters){
          kMeansRes = stats::kmeans(as.data.frame(mcoin$mcoa$SynVar), centers = i)[[1]]
          ClusterRes = data.frame(row.names = colnames(Data[[1]]),
                                Cluster = kMeansRes,
                                stringsAsFactors = FALSE)
          self$Multi_Omics$ClusterRes$MCIA[[paste0("Clusters_", i)]] = ClusterRes
          self$Multi_Omics$Fit$MCIA[[paste0("Clusters_", i)]] = mcoin
        }
      }else{
        self$Multi_Omics$Fit$MCIA[[paste0("Factors_", Components)]] = mcoin
      }

    },

    #jNMF
    #' @description
        #' Perform multi-omics integration using joint non-negative matrix
        #' factorization (jNMF) from the nnTensor package.
    #' @param Clusters The number of sample clusters to be created. Can be a
    #' single integer to calculate a single cluster assignment, or a range to
    #' calculate multiple cluster assignments. In jNMF, the number of clusters
    #' is the number of low dimensions or components created.
    #' @param NMFAlgorithm String indicating algorithm for the divergence
    #' between X and X-bar. Must be one of: `"KL"` (default), `"Frobenius"` or
    #' `"IS"`.
    #' @param MaxIter Number of iterations. Default is 100.
    #' @returns A list of results. Fit of the integration for each included
    #' cluster number, stored in $Multi_Omics$Fit$jNMF. Sample cluster
    #' assignments for each included cluster number, stored in
    #' $Multi_Omics$ClusterRes$jNMF. Component (factor) by feature matrices for
    #' each individual omics dataset, stored by included cluster number
    #' (= components in jNMF) in $Multi_Omics$FeatureRes$jNMF.
    #' @references Tsuyuzaki, K., and Nikaido, I. (2024). nnTensor: Non-Negative
    #' Tensor Decomposition. R package version 1.3.0
    run_jNMF = function(Clusters = NULL,
                        NMFAlgorithm = "KL", #algorithm for divergence between X and X_bar: frobenius, KL (default) or IS
                        MaxIter = 100 #number of iterations
    ){
      #data check
      #Same number of samples in columns needed for jNMF
      if(!is.null(self$Preprocessed_Omics$Matched_Data)){
        Data = self$Preprocessed_Omics$Matched_Data
      } else if(!is.null(self$Preprocessed_Omics$Data)){
        Data = self$Preprocessed_Omics$Data
      } else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      #non-negativity constraint, so add absolute value of smallest negative number to all values
      #also rescale each dataset to make magnitudes between sets comparable
      if(is.null(self$Preprocessed_Omics$Non_Neg_Data) == TRUE){
        self$Preprocessed_Omics$Non_Neg_Data = lapply(Data, function(dat){
          if(!all(dat >= 0)){
            dat = pmax(dat + abs(min(dat)), 0)
            #pmax is check, if any negatives remain, they should become 0
          }
          #scale data
          dat = dat/max(dat)
          #NMF methods require samples in rows, so transpose matrix
          dat = t(dat)
        })
      }

      Data = self$Preprocessed_Omics$Non_Neg_Data

      #run jNMF
      for(i in Clusters){
        jNMFres = nnTensor::jNMF(X = Data,
                                 pseudocount = .Machine$double.eps,
                                 J = i,
                                 algorithm = NMFAlgorithm,
                                 num.iter = MaxIter)
        #save cluster results
        jNMFclustres = max.col(jNMFres$W,
                                     ties.method = "first")

        ClusterRes = data.frame(Cluster = jNMFclustres,
                                row.names = rownames(Data[[1]]),
                                stringsAsFactors = FALSE)

        self$Multi_Omics$ClusterRes$jNMF[[paste0("Clusters_", i)]] = ClusterRes
        self$Multi_Omics$Fit$jNMF[[paste0("Clusters_", i)]] = jNMFres

        for(j in 1:length(Data)){
          jNMFFeatRes = as.data.frame(jNMFres$H[[j]])
          jNMFFeatRes = tibble::rownames_to_column(jNMFFeatRes, var = "Features")
          self$Multi_Omics$FeatureRes$jNMF[[paste0("Factors_", i)]][[j]] = jNMFFeatRes
        }
        names(self$Multi_Omics$FeatureRes$jNMF[[paste0("Factors_", i)]]) = names(Data)
      }
    },

    #iNMF
    #' @description
        #' Perform multi-omics integration and sample clustering using
        #' integrative non-negative matrix factorization from the IntNMF package.
    #' @param Clusters The number of sample clusters to be created. Can be a
    #' single integer to calculate a single cluster assignment, or a range to
    #' calculate multiple cluster assignments. In iNMF, the number of clusters
    #' is the number of low dimensions or components created.
    #' @param MaxIter Number of iterations. Default is 200.
    #' @param StabilityCount Count for stability in the connectivity matrix.
    #' Default is 20.
    #' @param MatrixInit Number of initializations of random matrices. Default
    #' is 30.
    #' @returns A list of results. Fit of the integration for each included
    #' cluster number, stored in $Multi_Omics$Fit$iNMF. Sample cluster
    #' assignments for each included cluster number, stored in
    #' $Multi_Omics$ClusterRes$iNMF. Component (factor) by feature matrices for
    #' each individual omics dataset, stroled by included cluster number
    #' (= components in iNMF) in $Multi_Omics$FeatureRes$iNMF.
    #' @references Chalise, P., Raghavan, R., and Fridley, B. (2025). IntNMF:
    #' Integrative Clustering of Multiple Genomic Dataset. R package version
    #' 1.3.0
    run_iNMF = function(Clusters = NULL,
                        MaxIter = 200,
                        StabilityCount = 20, #Count for stability in connectivity matrix
                        MatrixInit = 30 #number of random matrix initializations
    ){
      #data check
      ##Same number of samples in columns needed for iNMF
      if(!is.null(self$Preprocessed_Omics$Matched_Data)){
        Data = self$Preprocessed_Omics$Matched_Data
      } else if(!is.null(self$Preprocessed_Omics$Data)){
        Data = self$Preprocessed_Omics$Data
      } else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      #non-negativity constraint, so add absolute value of smallest negative number to all values
      #also rescale each dataset to make magnitudes between sets comparable
      if(is.null(self$Preprocessed_Omics$Non_Neg_Data) == TRUE){
        self$Preprocessed_Omics$Non_Neg_Data = lapply(Data, function(dat){
          if(!all(dat >= 0)){
            dat = pmax(dat + abs(min(dat)), 0)
          }
          #scale data
          dat = dat/max(dat)
          #NMF methods require samples in rows, so transpose matrix
          dat = t(dat)
        })
      }

      Data = self$Preprocessed_Omics$Non_Neg_Data

      for(i in Clusters){
        #run iNMF
        iNMFres = IntNMF::nmf.mnnals(dat = Data,
                                     k = i,
                                     maxiter = MaxIter,
                                     st.count = StabilityCount,
                                     n.ini = MatrixInit)
        iNMFclus = iNMFres$clusters

        #save clustering results
        ClusterRes = data.frame(Cluster = as.numeric(iNMFclus),
                                row.names = rownames(Data[[1]]),
                                stringsAsFactors = FALSE)

        self$Multi_Omics$ClusterRes$iNMF[[paste0("Clusters_", i)]] = ClusterRes
        self$Multi_Omics$Fit$iNMF[[paste0("Clusters_", i)]] = iNMFres

        for(j in 1:length(Data)){
          iNMFFeatRes = as.data.frame(t(iNMFres$H[[j]]))
          iNMFFeatRes = tibble::rownames_to_column(iNMFFeatRes, var = "Features")
          self$Multi_Omics$FeatureRes$iNMF[[paste0("Factors_", i)]][[j]] = iNMFFeatRes
        }
        names(self$Multi_Omics$FeatureRes$iNMF[[paste0("Factors_", i)]]) = names(Data)
      }
    },

    #LRAcluster
    #' @description
        #' Perform multi-omics integration and sample clustering using
        #' LRAcluster for R as done in the MOVICS package.
    #' @param Clusters The Number of sample clusters to be created. Can be a
    #' single integer to calculate a single cluster assignment, or a range to
    #' calculate multiple cluster assignments.
    #' @param Linkage Agglomeration method to be used for the hierarchical
    #' clustering. Must be one of: `"ward.D"` (default), `"ward.D2"`, `"single"`,
    #' `"complete"`, `"average"`, `"mcquitty"`, `"median"` or `"centroid"`.
    #' @returns A list of results. Fit of the integration for each included
    #' cluster number, stored in $Multi_Omics$Fit$LRAcluster. Sample cluster
    #' assignments for each included cluster number, stored in
    #' $Multi_Omics$ClusterRes$LRAcluster. Matrix with coordinates of all
    #' samples in reduced space for each included cluster number, stored in
    #' $Multi_Omics$CoordData$LRAcluster.
    #' @references Wu D, Wang D, Zhang MQ, Gu J (2015). Fast dimension reduction
    #' and integrative clustering of multi-omics data using low-rank approximation:
    #' application to cancer molecular classification. BMC Genomics, 16(1):1022.
    #' @references   * Lu, X., Meng, J., Zhou, Y., Jiang, L., and Yan, F. (2020).
    #'   MOVICS: an R package for multi-omics integration and visualization in
    #'   cancer subtyping. Bioinformatics, btaa1018.
    #'   * Wu D, Wang D, Zhang MQ, Gu J (2015). Fast dimension reduction and
    #'   integrative clustering of multi-omics data using low-rank approximation:
    #'   application to cancer molecular classification. BMC Genomics, 16(1):1022.
    run_LRAcluster = function(Clusters = NULL,
                              Linkage = "ward.D" #hclust agglomeration method: ward.D,ward.D2, single, complete (default), average, mcquitty, median, centroid
    ){

      #data check
      ##Same number of samples in columns needed for LRAcluster
      if(!is.null(self$Preprocessed_Omics$Matched_Data)){
        Data = self$Preprocessed_Omics$Matched_Data
      } else if(!is.null(self$Preprocessed_Omics$Data)){
        Data = self$Preprocessed_Omics$Data
      } else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      for(i in Clusters){
        #run getLRAcluster
        LRAres = getLRAcluster_MOVICS(data = Data,
                                      N.clust = i,
                                      type = rep("gaussian", length(Data)),
                                      clusterAlg = Linkage)

        #save cluster results
        LRAres$clust.res = subset(LRAres$clust.res, select = -samID)
        names(LRAres$clust.res)[names(LRAres$clust.res) == "clust"] <- "Cluster"
        ClusterRes = LRAres$clust.res

        self$Multi_Omics$ClusterRes$LRAcluster[[paste0("Clusters_", i)]] = ClusterRes
        self$Multi_Omics$Fit$LRAcluster[[paste0("Clusters_", i)]] = LRAres$fit
        Coord = as.data.frame(t(LRAres$fit$coordinate))
        rownames(Coord) = colnames(Data[[1]])
        self$Multi_Omics$CoordData$LRAcluster[[paste0("Factors_", i)]] = Coord
      }
    },

    #COCA
    #' @description
        #' Perform multi-omics integration and sample clustering using cluster
        #' of cluster analysis from the COCA package.
    #' @param DatasetClusters Number of clusters per omics dataset. If a single
    #' integer, all omics datasets get the same number of clusters. If a vector
    #' with the same length as the number of omics datasets, each dataset can get
    #' a different number of clusters. If NULL, the number of clusters per
    #' dataset will be estimated using the silhouette method, and argument
    #' `"maxClusters` should be included.
    #' @param maxClusters Maximum number of clusters per omics dataset if
    #' `DatasetClusters` is NULL. If a single integer, all omics datasets have the
    #' same maximum number of clusters considered. If a vector with the same
    #' length as the number of omics datasets, each dataset has an individual
    #' maximum number of clusters considered. Default is 12.
    #' @param DatasetClusterMethod String indicating the method to cluster the
    #' individual omics datasets. If single string, all omics datasets are
    #' clustered with the same method. If vector of strings, each omics dataset
    #' is clustered with an individual method. Must be one of: `"kmeans"`,
    #' `"hclust"` (default) or `"pam"`.
    #' @param DatasetDistance Distance to be used in clustering of the
    #' individual omics datasets. If single string, all omics datasets are
    #' clustered with the same distance. If vector of strings, each omics
    #' dataset is clustered with an individual distances. Must be one of:
    #' `"euclidean"` (default, available for `"kmeans"`, `"hclust"` and `"pam"`),
    #' `"manhattan"` (available for `"kmeans"`, `"hclust"` and `"pam"`),
    #' `"gower"` (available for `"pam"`), `"maximum"` (available for `"kmeans"`
    #' and `"hclust"`), `"canberra"` (available for `"kmeans"` and `"hclust"`),
    #' `"binary"` (available for `"kmeans"` and `"hclust"`) or `"minkowski"`
    #' (available for `"kmeans"` and `"hclust"`).
    #' @param FullClusters Number of clusters after integration. Can be a single
    #' integer to calculate a single cluster assignment, or a range to calculate
    #' multiple cluster assignments. If NULL, the number of clusters will be
    #' estimated using the silhouette method, and argument `"maxFullClusters"`
    #' should be included.
    #' @param maxFullClusters Maximum number of clusters if `FullClusters` is
    #' NULL. Default is 12.
    #' @param FullHClustMethod Agglomeration method to be used for the
    #' hierarchical clustering of the integrated dataset. Must be one of:
    #' `"ward.D"`, `"ward.D2"`, `"single"`, "complete"`, `"average"` (default),
    #' `"mcquitty"`, `"median"` or`"centroid"`.
    #' @param FullClusterMethod String indicating the method to cluster the
    #' integrated dataset. Must be one of: `"kmeans"` (default) or `"hclust"`.
    #' @param FullDistance Distance to be used in clustering of the
    #' integrated dataset. Must be one of: `"pearson"`, `"spearman"`,
    #' `"euclidean"` (default), `"manhattan"`, `"maximum"`, `"canberra"`,
    #' `"binary"` or `"minkowski"`.
    #' @returns A list of results. Fit of the integration for each included
    #' `FullClusters` cluster number, stored in $Multi_Omics$Fit$COCA. Sample
    #' cluster assignments for each included `FullClusters` cluster number,
    #' stored in $Multi_Omics$ClusterRes$COCA.
    #' @references Cabassi A, Kirk PD (2020). “Multiple kernel learning for
    #' integrative consensus clustering of genomic datasets.” Bioinformatics.
    run_COCA =  function(DatasetClusters = NULL, #number of clusters per dataset (integer ==> all datasets have same number of clusters, vector ==> per dataset)
                         maxClusters = 12, #maximum number of clusters to be considered per dataset
                         DatasetClusterMethod = "hclust", #clustering method, vector length 1 ==> used for all datasets, vector length(Data) ==> specific method per dataset. k-means, hclust (default) or pam (partitioning around medoids)
                         DatasetDistance = "euclidean", #distances used for clustering steps,vector length 1 ==> used for all datasets, vector length(Data) ==> specific method per dataset. euclidean (default), manhattan, gower (pam only), maximum (not pam), canberra (not                           pam), binary (not pam), minkowski (not pam)
                         FullClusters = NULL, #number of clusters for cluster of clusters analysis
                         maxFullClusters = 12, #maximum number of clusters to be considered for cluster of clusters analysis
                         FullHClustMethod = "average", #hclust agglomeration method: ward.D, ward.D2, single, complete (default), average, mcquitty, median, centroid
                         FullClusterMethod = "kmeans", #clustering inside of consensus clustering, kmeans (default) or hclust
                         FullDistance = "euclidean" #distance for CC hclust, can be pearson, spearman, euclidean (default), maximum, manhattan, canberra, binary, minkowski
    ){

      #data checks
      #No NA, so impute missings if needed
      if(!is.null(self$Preprocessed_Omics$Matched_Data)){
        Data = self$Preprocessed_Omics$Matched_Data
      } else if(!is.null(self$Preprocessed_Omics$Data)){
        Data = self$Preprocessed_Omics$Data
      } else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      #samples required in rows
      tData = lapply(Data, t)

      #parameter for calculating optimal cluster number for dataset calculations
      if(is.null(DatasetClusters) == TRUE){
        DatasetOptCluster = TRUE
      } else{
        DatasetOptCluster = FALSE
      }

      #parameter for calculating optimal cluster number for cluster of clusters analysis
      if(is.null(FullClusters) == TRUE){
        FullOptCluster = TRUE
      } else{
        FullOptCluster = FALSE
      }

      #build matrix of clusters
      outputMOC = coca::buildMOC(data = tData,
                                 M = length(tData), #number of datasets
                                 K = DatasetClusters, #if NULL, cluster number will be calculated for each dataset based on maxK
                                 maxK = maxClusters,
                                 methods = DatasetClusterMethod,
                                 distances = DatasetDistance,
                                 fill = FALSE,
                                 computeAccuracy = FALSE,
                                 fullData = FALSE,
                                 widestGap = DatasetOptCluster,
                                 dunns = DatasetOptCluster,
                                 dunn2s = DatasetOptCluster
      )
      moc = outputMOC$moc

      if(is.null(FullClusters)){
        coca = coca::coca(moc = moc,
                          K = FullClusters, #if NULL, cluster number will be calculated based on maxK
                          maxK = maxFullClusters,
                          hclustMethod = FullHClustMethod,
                          choiceKmethod = "silhouette",
                          ccClMethod = FullClusterMethod,
                          ccDistHC = FullDistance,
                          widestGap = FullOptCluster,
                          dunns = FullOptCluster,
                          dunn2s = FullOptCluster
        )

        #save clustering results
        ClusterRes = data.frame(row.names = colnames(Data[[1]]),
                                Cluster = as.numeric(coca$clusterLabels),
                                stringsAsFactors = FALSE)
        NClust = max(as.numeric(coca$clusterLabels))

        self$Multi_Omics$ClusterRes$COCA[[paste0("Clusters_", NClust)]] = ClusterRes
        self$Multi_Omics$Fit$COCA[[paste0("Clusters_", NClust)]] = c(outputMOC, coca)
      }else{
        for(i in FullClusters){
          coca = coca::coca(moc = moc,
                            K = i, #if NULL, cluster number will be calculated based on maxK
                            maxK = maxFullClusters,
                            hclustMethod = FullHClustMethod,
                            choiceKmethod = "silhouette",
                            ccClMethod = FullClusterMethod,
                            ccDistHC = FullDistance,
                            widestGap = FullOptCluster,
                            dunns = FullOptCluster,
                            dunn2s = FullOptCluster
          )

          #save clustering results
          ClusterRes = data.frame(row.names = colnames(Data[[1]]),
                                  Cluster = as.numeric(coca$clusterLabels),
                                  stringsAsFactors = FALSE)
          NClust = max(as.numeric(coca$clusterLabels))

          self$Multi_Omics$ClusterRes$COCA[[paste0("Clusters_", NClust)]] = ClusterRes
          self$Multi_Omics$Fit$COCA[[paste0("Clusters_", NClust)]] = c(outputMOC, coca)
        }
      }

    },

    #MOFA
    #' @description
        #' Perform multi-omics integration using multi-omics factor analysis
        #' from the MOFA2 package.
    #' @param PerformCluster Whether to perform k-means clustering on the sample
    #' embeddings. Default is `TRUE`.
    #' @param Clusters The number of sample clusters to be created. Can be a
    #' single integer to calculate a single cluster assignment, or a range to
    #' calculate multiple cluster assignments.
    #' @param Factors Number of factors (components).
    #' @param ScaleViews Whether to scale the different omics datasets to be
    #' scaled to have the same unit variance. Default is `FALSE`.
    #' @param SpikeSlabFactors Whether to use spike and slab sparsity on the
    #' factors (components). Default is `FALSE`.
    #' @param SpikeSlabWeights Whether to use spike and slab sparsity on the
    #' features. Default is `TRUE`.
    #' @param ARDFactors Whether to use ARD sparsity on the factors. Default is
    #' `FALSE`.
    #' @param ARDWeights Whether to use ARD sparsity on the weights. Default is
    #' `TRUE`.
    #' @param MaxIter Maximum number of iterations. Default is 1000.
    #' @param Convergence String indicating the convergence criteria. Must be
    #' one of: `"fast"` (default), `"medium"`, or `"slow"`.
    #' @param StartELBO Integer indicating the first iteration to compute the
    #' ELBO statistic. Default is 1.
    #' @param FreqELBO Integer indicating the frequency of ELBO computations.
    #' Default is 1.
    #' @param Stochastic Whether to use stochastic inference. Default is `FALSE`.
    #' @param SaveData Whether to save the trained MOFA model. Default is
    #' `FALSE`.
    #' @param Outfile File path to save the trained MOFA model.
    #' @returns A list of results. Fit of the integration model, both pretrained
    #' and trained, stored by number of factors in $Multi_Omics$Fit$MOFA. Factor
    #' (component) by feature matrices with feature weights for each omics
    #' dataset, stored by number of factors in $Multi_Omics$FeatureRes$MOFA.
    #' Sample by factor (component) matrix, stored by number of factors in
    #' $Multi_Omics$FactorRes$MOFA. Total variance explained for each factor
    #' (component) per omics type, stored by number of factors in
    #' $Multi_Omics$VarianceExp$MOFA. Matrix with coordinates of all
    #' samples in reduced factor (component) space, stored by number of factors
    #' in $Multi_Omics$CoordData$MOFA.
    #' @references Argelaguet R, Velten B, Arnol D, Dietrich S, Zenz T, Marioni
    #' JC, Buettner F, Huber W, Stegle O (2018). “Multi‐Omics Factor Analysis—a
    #' framework for unsupervised integration of multi‐omics data sets.”
    #' Molecular Systems Biology, 14.
    run_MOFA = function(PerformCluster = TRUE, #whether to perform clustering
                        Clusters = NULL, #number of clusters
                        Factors = NULL,
                        ScaleViews = FALSE, #scale each view to unit variance
                        SpikeSlabFactors = FALSE, #spike-slab sparsity prior in factors
                        SpikeSlabWeights = TRUE, #spike-slab sparsity prior in weights
                        ARDFactors = FALSE, #ARD prior in factors
                        ARDWeights = TRUE, #ARD prior in weights
                        MaxIter = 1000,
                        Convergence = "fast", #options: "fast" (default), "medium", "slow" (fast for exploration)
                        StartELBO = 1, #initial iteration to compute ELBO
                        FreqELBO = 1, #frequency of ELBO computation
                        Stochastic = FALSE, #stochastic interference
                        SaveData = FALSE, #whether to save trained MOFA model
                        Outfile = NULL #where to save trained MOFA model
    ){
      #data checks
      ##Same number of samples in columns needed for MOFA
      if(!is.null(self$Preprocessed_Omics$Matched_Data)){
        Data = self$Preprocessed_Omics$Matched_Data
      } else if(!is.null(self$Preprocessed_Omics$Data)){
        Data = self$Preprocessed_Omics$Data
      } else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      #create MOFA object
      MOFAobject = MOFA2::create_mofa(data = Data,
                                      groups = NULL)

      #prepare MOFA
      data_opts = MOFA2::get_default_data_options(MOFAobject)
      model_opts = MOFA2::get_default_model_options(MOFAobject)
      train_opts = MOFA2::get_default_training_options(MOFAobject)

      #set data options
      data_opts$scale_views = ScaleViews
      data_opts$scale_groups = FALSE
      data_opts$center_groups = TRUE

      #set model options
      model_opts$likelihoods = rep("gaussian", length(Data))
      model_opts$num_factors = Factors
      model_opts$spikeslab_factors = SpikeSlabFactors
      model_opts$spikeslab_weights = SpikeSlabWeights
      model_opts$ard_factors = ARDFactors
      model_opts$ard_weights = ARDWeights

      #set training options
      train_opts$maxiter = MaxIter
      train_opts$convergence_mode = Convergence
      train_opts$startELBO = StartELBO
      train_opts$freqELBO = FreqELBO
      train_opts$stochastic = Stochastic

      MOFAobject = MOFA2::prepare_mofa(object = MOFAobject,
                                       data_options = data_opts,
                                       model_options = model_opts,
                                       training_options = train_opts)

      #train MOFA model
      MOFAmodel = MOFA2::run_mofa(object = MOFAobject,
                                  outfile = Outfile,
                                  use_basilisk = TRUE,
                                  save_data = SaveData)

      #save feature data
      MOFAFeatRes = lapply(MOFAmodel@expectations$W, as.data.frame)
      for(i in 1:length(Data)){
        MOFAFeats = tibble::rownames_to_column(MOFAFeatRes[[i]], var = "Features")
        self$Multi_Omics$FeatureRes$MOFA[[paste0("Factors_", Factors)]][[i]] = MOFAFeats
      }
      names(self$Multi_Omics$FeatureRes$MOFA[[paste0("Factors_", Factors)]]) = names(Data)
      #save factor data
      self$Multi_Omics$FactorRes$MOFA[[paste0("Factors_", Factors)]] = lapply(MOFAmodel@expectations$Z, as.data.frame)

      self$Multi_Omics$VarianceExp$MOFA[[paste0("Factors_", Factors)]] = MOFAmodel@cache$variance_explained
      Coord = as.data.frame(MOFAmodel@expectations$Z$group1)
      rownames(Coord) = colnames(Data[[1]])
      self$Multi_Omics$CoordData$MOFA[[paste0("Factors_", Factors)]] = Coord


      #create clusters and save results
      if(PerformCluster == TRUE){
        for(i in Clusters){
          kMeansRes = stats::kmeans(as.data.frame(MOFAmodel@expectations$Z$group1), centers = i)[[1]]
          ClusterRes = data.frame(row.names = colnames(Data[[1]]),
                                  Cluster = kMeansRes,
                                  stringsAsFactors = FALSE)
          self$Multi_Omics$ClusterRes$MOFA[[paste0("Clusters_", i)]] = ClusterRes
          self$Multi_Omics$Fit$MOFA$Pretrained[[paste0("Clusters_", i)]] = MOFAobject
          self$Multi_Omics$Fit$MOFA$Model[[paste0("Clusters_", i)]] = MOFAmodel
        }
      }else{
        self$Multi_Omics$Fit$MOFA$Pretrained[[paste0("Factors_", Factors)]] = MOFAobject
        self$Multi_Omics$Fit$MOFA$Model[[paste0("Factors_", Factors)]] = MOFAmodel
        self$Multi_Omics$VarianceExp$MOFA[[paste0("Factors_", Factors)]] = MOFAmodel@cache$variance_explained
      }


    },

    #SNF
    #' @description
        #' Perform multi-omics integration and sample clustering using
        #' similarity network fusion from the SNFtool package.
    #' @param Distance String indicating the distance used to calculate the
    #' distance matrices for each omics dataset. If single string, distance
    #' matrices are only calculated for that distance. If a vector of strings,
    #' distance matrices are calculated for all distances, and clustering will
    #' be performed for all distances. Can be: `"euclidean squared"` (default
    #' in SNF), `"euclidean"`, `"manhattan"`, `"minkowski 0.25"`,
    #' `"minkowski 0.5"`, `"minkowski 3"`, `"minkowski 4"`.
    #' @param Clusters Number of sample clusters to be created. Can be a single
    #' integer to calculate a single cluster assignment, or a range to calculate
    #' multiple cluster assignments.
    #' @param Neighbors Number of nearest neighbors. Default is 20.
    #' @param Sigma Variance for the local model. Default is 0.5.
    #' @param Iter Number of iterations for the diffusion process. Default is 10.
    #' @param SpecClusterType Variant of spectral clustering to use. Must be one
    #' of: 1, 2 or 3 (default).
    #' @returns A list of results. Fit of the integration for each included
    #' cluster number and distance metric, stored in $Multi_Omics$Fit$SNF.
    #' Sample cluster assignments for each included cluster number and distance
    #' metric, stored in $Multi_Omics$ClusterRes$SNF.
    #' @references Wang, B., Mezlini, A., Demir, F., Fiume, M., Tu, Z., Brudno,
    #' M., Haibe-Kains, B., and Goldenberg, A. (2021) SNFtool: Similarity
    #' Network Fusion. R package version 2.3.1
    run_SNF = function(Distance = "euclidean squared",
                       Clusters = NULL,
                       Neighbors = 20, #number of neighbours for affinity matrix (between 10-30)
                       Sigma = 0.5, #hyperparameter for affinity matrix (between 0.3-0.8)
                       Iter = 10, #number of iterations (between 10-20)
                       SpecClusterType = 3 #type of spectral clustering (1-3, default 3)
    ){
      #data check
      if(!is.null(self$Preprocessed_Omics$Matched_Data)){
        Data = self$Preprocessed_Omics$Matched_Data
      } else if(!is.null(self$Preprocessed_Omics$Data)){
        Data = self$Preprocessed_Omics$Data
      } else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      #requires samples in rows, so transpose Data
      Data = lapply(Data, t)

      #run SNF
      Distance_List = Distance
      for(distance in Distance_List){
        AffMatrix = SNF_Distance_Affinity_Matrix(Distance = distance,
                                                 Data = Data,
                                                 Neighbours = Neighbors,
                                                 Sigma = Sigma)
        W = SNFtool::SNF(Wall = AffMatrix,
                         K = Neighbors,
                         t = Iter)
        for(i in Clusters){
          clust = SNFtool::spectralClustering(W,
                                              i,
                                              type = SpecClusterType)
          ClusterRes = data.frame(row.names = rownames(Data[[1]]),
                                  Cluster = clust,
                                  stringsAsFactors = FALSE)
          self$Multi_Omics$Fit$SNF[[distance]][[paste0("Clusters_", i)]] = W
          self$Multi_Omics$ClusterRes$SNF[[distance]][[paste0("Clusters_", i)]] = ClusterRes
        }
      }
    },

    #' @description
        #' Multi-omics integration and sample clustering using UMAP dimensionality
        #' reduction as described in the GAUDI method. In addition to clustering
        #' with HDBSCAN, as done in GAUDI, we provide DBSCAN and k-means
        #' methods for sample clustering.
    #' @param SingleOmicsUMAPNeighbors Size of the local neighborhood used for
    #' manifold approximation for each individual omics dataset. If a single
    #' integer, all omics datasets get the same size. If a vector with the same
    #' length as the number of omics datasets, each omics dataset gets a
    #' different size. Default is 15, recommended range is 2-100.
    #' @param SingleOmicsUMAPDimensions The number of dimensions into which the
    #' individual omics datasets are embedded. If a single integer, all omics
    #' datasets get the same number of dimensions. If a vector with the same
    #' length as the number of omics datasets, each omics dataset gets a
    #' different number of dimensions. Default is 4.
    #' @param SingleOmicsUMAPDistance The distance metric used to find the
    #' nearest neighbors in the individual omics datasets. If a single string,
    #' all omics datasets use the same distance metric. If a vector of strings
    #' with the same length as the number of omics datasets, each omics dataset
    #' uses a different distance metric. Must be one of: `"euclidean"`
    #' (default), `"cosine"`, `"manhattan"`, `"hamming"`, `"correlation"` or
    #' `"categorical"`.
    #' @param SingleOmicsUMAPScale Whether to perform scaling on the features of
    #' the individual omics datasets. If a single boolean, all omics datasets
    #' use the same scaling setting. If a vector of booleans, each omics dataset
    #' uses a different scaling setting. Default is `FALSE`.
    #' @param SingleOmicsUMAPInit Type of initialization for the coordinates of
    #' the individual omics datasets. If a single string, all omics datasets use
    #' the same initialization type. If a vector of strings, each omics dataset
    #' uses a different initialization type. Must be one of: `"spectral"`
    #' (default), `"normlaplacian"`, `"random"`, `"lvrandom"`, `"laplacian"`,
    #' `"pca"`, `"spca"` or `"agspectral"`.
    #' @param SingleOmicsUMAPPCA Whether to reduce the number of features in the
    #' individual omics datasets using PCA to increase the performance. May come
    #' at a cost to accuracy. A positive integer less than the number of
    #' features in the omics dataste is required. If a single integer, all omics
    #' datasets are reduced to the same number of principle components. If a
    #' vector of integers, each omics dataset is reduced to a different number
    #' of principle components. Default is `NULL`.
    #' @param ConcatUMAPNeighbors Size of the local neighborhood used for
    #' manifold approximation for the UMAP calculation of the concatenated
    #' single omics UMAP embeddings. Default is 15, recommended range is 2-100.
    #' @param ConcatUMAPDimensions The number of dimensions into which the
    #' concatenated single omics UMAP embeddings are embedded. Default is 4.
    #' @param ConcatUMAPDistance The distance metric used to find the nearest
    #' neighbors for the concatenated single omics UMAP embeddings. Must be one
    #' of: `"euclidean"` (default), `"cosine"`, `"manhattan"`, `"hamming"`,
    #' `"correlation"` or `"categorical"`.
    #' @param ConcatUMAPScale Whether to perform scaling on the UMAP dimensions
    #' of the concatenated single omics UMAP embeddings. Default is `FALSE`.
    #' @param ConcatUMAPInit Type of initialization for the coordinates of the
    #' concatenated single omics UMAP embeddings. Must be one of: `"spectral"`
    #' (default), `"nornlaplacian"`, `"random"`, `"lvrandom"`, `"laplacian"`,
    #' `"pca"`, `"spca"` or `"agspectral"`.
    #' @param ConcatUMAPPCA Whether to reduce the number of UMAP dimensions of
    #' the concatenated single omics UMAP embeddings using PCA to increase the
    #' performance. May come at a cost to accuracy. Requires a positive integer
    #' less than the the total number of UMAP dimension in all single omics UMAP
    #' embeddings. Default is `NULL`.
    #' @param ClusterMethod Method used for clustering of the concatenated UMAP
    #' embedding. Must be one of: `"HDBSCAN"`, `"DBSCAN"` or `"Kmeans"`.
    #' @param HDBSCANminPts Integer of the minimal size of the clusters. Default
    #' is 2. Only included if `ClusterMethod` = `"HDBSCAN"`.
    #' @param DBSCANeps Radius of the epsilon neighborhood. Default is `NULL`.
    #' Only included if `ClusterMethod` = `"DBSCAN"`.
    #' @param DBSCANminPts Integer of the minimal size of the clusters. Default
    #' is 2. Only included if `ClusterMethod` = `"DBSCAN"`.
    #' @param KmeansClusters Number of sample clusters to be created. Can be a
    #' single integer to calculate a single cluster assignment, or a range to
    #' calculate multiple cluster assignments. Only included if `ClusterMethod`
    #' = `"Kmeans"`.
    #' @returns A list of results. Fit of the single omics UMAPs, concatenated
    #' UMAP and clustering method for each cluster number, stored by clustering
    #' method in $Multi_Omics$Fit$GAUDI. Sample cluster assignments for each
    #' cluster number, stored by clustering method in
    #' $Multi_Omics$ClusterRes$GAUDI. Sample coordinates in the concatenated
    #' UMAP embedding for each cluster number, stored by clustering method in
    #' $Multi_Omics$CoordData$GAUDI.
    #' @references   * Castellano-Escuder P, Zachman DK, Han K, Hirschey MD.
    #'   GAUDI: interpretable multi-omics integration with UMAP embeddings and
    #'   density-based clustering. Nat Commun. 2025 Jul 1;16(1):5771.
    #'   * Hahsler M, Piekenbrock M (2025). dbscan: Density-Based Spatial
    #'    Clustering of Applications with Noise (DBSCAN) and Related Algorithms.
    #'    doi:10.32614/CRAN.package.dbscan, R package version 1.2.3
    run_GAUDI = function(SingleOmicsUMAPNeighbors = 15,
                         SingleOmicsUMAPDimensions = 4,
                         SingleOmicsUMAPDistance = "euclidean",
                         SingleOmicsUMAPScale = FALSE,
                         SingleOmicsUMAPInit = "spectral",
                         SingleOmicsUMAPPCA = NULL,
                         ConcatUMAPNeighbors = 15,
                         ConcatUMAPDimensions = 4,
                         ConcatUMAPDistance = "euclidean",
                         ConcatUMAPScale = FALSE,
                         ConcatUMAPInit = "spectral",
                         ConcatUMAPPCA = NULL,
                         ClusterMethod = "HDBSCAN", #method for clustering, options: HDBSCAN, DBSCAN, Kmeans
                         HDBSCANminPts = 2, #minimum cluster size for HDBSCAN
                         DBSCANeps = NULL, #radius of epsilon neighborhood
                         DBSCANminPts = 2, #minimum cluster size for DBSCAN
                         KmeansClusters = NULL #number of clusters for K-means
    ){
      #store parameters
      CurrentParam = c(SingleOmicsUMAPNeighbors = SingleOmicsUMAPNeighbors, SingleOmicsUMAPDimensions = SingleOmicsUMAPDimensions, SingleOmicsUMAPDistance = SingleOmicsUMAPDistance, SingleOmicsUMAPScale = SingleOmicsUMAPScale, SingleOmicsUMAPInit = SingleOmicsUMAPInit, SingleOmicsUMAPPCA = SingleOmicsUMAPPCA, ConcatUMAPNeighbors = ConcatUMAPNeighbors, ConcatUMAPDimensions = ConcatUMAPDimensions, ConcatUMAPDistance = ConcatUMAPDistance, ConcatUMAPScale = ConcatUMAPScale, ConcatUMAPInit = ConcatUMAPInit, ConcatUMAPPCA = ConcatUMAPPCA)

      #Data selection
      if(!is.null(self$Preprocessed_Omics$Matched_Data)){
        Data = self$Preprocessed_Omics$Matched_Data
      } else if(!is.null(self$Preprocessed_Omics$Data)){
        Data = self$Preprocessed_Omics$Data
      } else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      #UMAP requires samples in rows, so transpose data
      tData = lapply(Data, t)

      #UMAP
      if(!is.null(self$Multi_Omics$Fit$GAUDI)){
        GAUDIClustMethod = names(self$Multi_Omics$Fit$GAUDI)
        for(method in GAUDIClustMethod){
          GAUDIClustNum = names(self$Multi_Omics$Fit$GAUDI[[method]])
          for(clust in GAUDIClustNum){
            if(identical(self$Multi_Omics$Fit$GAUDI[[method]][[clust]]$Param, CurrentParam)){
              Single_Omics_UMAPs = self$Multi_Omics$Fit$GAUDI[[method]][[clust]]$Single_Omics_UMAPs
              Concat_UMAP = self$Multi_Omics$Fit$GAUDI[[method]][[clust]]$Concatenated_UMAP
              break
            }else{
              Single_Omics_UMAPs = list()

              #single omics UMAP
              #run each omics dataset with own parameters
              if(length(SingleOmicsUMAPNeighbors) | length(SingleOmicsUMAPDimensions) | length(SingleOmicsUMAPDistance) | length(SingleOmicsUMAPScale) | length(SingleOmicsUMAPInit) | length(SingleOmicsUMAPPCA) > 1){
                if(length(SingleOmicsUMAPNeighbors) == 1){
                  SingleOmicsUMAPNeighbors = rep(SingleOmicsUMAPNeighbors, each = length(tData))
                }
                if(length(SingleOmicsUMAPDimensions) == 1){
                  SingleOmicsUMAPDimensions = rep(SingleOmicsUMAPDimensions, each = length(tData))
                }
                if(length(SingleOmicsUMAPDistance) == 1){
                  SingleOmicsUMAPDistance = rep(SingleOmicsUMAPDistance, each = length(tData))
                }
                if(length(SingleOmicsUMAPScale) == 1){
                  SingleOmicsUMAPScale = rep(SingleOmicsUMAPScale, each = length(tData))
                }
                if(length(SingleOmicsUMAPInit) == 1){
                  SingleOmicsUMAPInit = rep(SingleOmicsUMAPInit, each = length(tData))
                }
                if(length(SingleOmicsUMAPPCA) == 1){
                  SingleOmicsUMAPPCA = rep(SingleOmicsUMAPPCA, each = length(tData))
                }

                Omics_Names = names(tData)
                for(i in 1:length(tData)){
                  Single_Omics_UMAPs[[i]] = uwot::umap(tData[[i]],
                                                       n_neighbors = SingleOmicsUMAPNeighbors[i],
                                                       n_components = SingleOmicsUMAPDimensions[i],
                                                       metric = SingleOmicsUMAPDistance[i],
                                                       scale = SingleOmicsUMAPScale[i],
                                                       init = SingleOmicsUMAPInit[i],
                                                       pca = SingleOmicsUMAPPCA[i])
                }
                names(Single_Omics_UMAPs) = Omics_Names
              } else{
                Single_Omics_UMAPs = lapply(tData,
                                            uwot::umap,
                                            n_neighbors = SingleOmicsUMAPNeighbors,
                                            n_components = SingleOmicsUMAPDimensions,
                                            metric = SingleOmicsUMAPDistance,
                                            scale = SingleOmicsUMAPScale,
                                            init = SingleOmicsUMAPInit,
                                            pca = SingleOmicsUMAPPCA)
              }

              for(i in length(Single_Omics_UMAPs)){
                colnames(Single_Omics_UMAPs[[i]])[1:ncol(Single_Omics_UMAPs[[i]])] = paste0("UMAP", 1:ncol(Single_Omics_UMAPs[[i]]))
              }

              #Single Omics UMAP concatenation
              SO_UMAP_Concat = dplyr::bind_cols(Single_Omics_UMAPs, .name_repair = "unique_quiet")

              #Concatenated data UMAP
              Concat_UMAP = uwot::umap(SO_UMAP_Concat,
                                       n_neighbors = ConcatUMAPNeighbors,
                                       n_components = ConcatUMAPDimensions,
                                       metric = ConcatUMAPDistance,
                                       scale = ConcatUMAPScale,
                                       init = ConcatUMAPInit,
                                       pca = ConcatUMAPPCA)

              colnames(Concat_UMAP)[1:ncol(Concat_UMAP)] = paste0("UMAP", 1:ncol(Concat_UMAP))
            }
          }
        }
      }else{
        Single_Omics_UMAPs = list()

        #single omics UMAP
        #run each omics dataset with own parameters
        if(length(SingleOmicsUMAPNeighbors) | length(SingleOmicsUMAPDimensions) | length(SingleOmicsUMAPDistance) | length(SingleOmicsUMAPScale) | length(SingleOmicsUMAPInit) | length(SingleOmicsUMAPPCA) > 1){
          if(length(SingleOmicsUMAPNeighbors) == 1){
            SingleOmicsUMAPNeighbors = rep(SingleOmicsUMAPNeighbors, each = length(tData))
          }
          if(length(SingleOmicsUMAPDimensions) == 1){
            SingleOmicsUMAPDimensions = rep(SingleOmicsUMAPDimensions, each = length(tData))
          }
          if(length(SingleOmicsUMAPDistance) == 1){
            SingleOmicsUMAPDistance = rep(SingleOmicsUMAPDistance, each = length(tData))
          }
          if(length(SingleOmicsUMAPScale) == 1){
            SingleOmicsUMAPScale = rep(SingleOmicsUMAPScale, each = length(tData))
          }
          if(length(SingleOmicsUMAPInit) == 1){
            SingleOmicsUMAPInit = rep(SingleOmicsUMAPInit, each = length(tData))
          }
          if(length(SingleOmicsUMAPPCA) == 1){
            SingleOmicsUMAPPCA = rep(SingleOmicsUMAPPCA, each = length(tData))
          }

          Omics_Names = names(tData)
          for(i in 1:length(tData)){
            Single_Omics_UMAPs[[i]] = uwot::umap(tData[[i]],
                                                 n_neighbors = SingleOmicsUMAPNeighbors[i],
                                                 n_components = SingleOmicsUMAPDimensions[i],
                                                 metric = SingleOmicsUMAPDistance[i],
                                                 scale = SingleOmicsUMAPScale[i],
                                                 init = SingleOmicsUMAPInit[i],
                                                 pca = SingleOmicsUMAPPCA[i])
          }
          names(Single_Omics_UMAPs) = Omics_Names
        } else{
          Single_Omics_UMAPs = lapply(tData,
                                      uwot::umap,
                                      n_neighbors = SingleOmicsUMAPNeighbors,
                                      n_components = SingleOmicsUMAPDimensions,
                                      metric = SingleOmicsUMAPDistance,
                                      scale = SingleOmicsUMAPScale,
                                      init = SingleOmicsUMAPInit,
                                      pca = SingleOmicsUMAPPCA)
        }

        for(i in length(Single_Omics_UMAPs)){
          colnames(Single_Omics_UMAPs[[i]])[1:ncol(Single_Omics_UMAPs[[i]])] = paste0("UMAP", 1:ncol(Single_Omics_UMAPs[[i]]))
        }

        #Single Omics UMAP concatenation
        SO_UMAP_Concat = dplyr::bind_cols(Single_Omics_UMAPs, .name_repair = "unique_quiet")

        #Concatenated data UMAP
        Concat_UMAP = uwot::umap(SO_UMAP_Concat,
                                 n_neighbors = ConcatUMAPNeighbors,
                                 n_components = ConcatUMAPDimensions,
                                 metric = ConcatUMAPDistance,
                                 scale = ConcatUMAPScale,
                                 init = ConcatUMAPInit,
                                 pca = ConcatUMAPPCA)

        colnames(Concat_UMAP)[1:ncol(Concat_UMAP)] = paste0("UMAP", 1:ncol(Concat_UMAP))
      }

      #clustering
      #HDBSCAN
      if(ClusterMethod == "HDBSCAN"){
        UMAP_clustres = dbscan::hdbscan(Concat_UMAP, minPts = HDBSCANminPts)[[1]]

        #save results
        ClusterRes = data.frame(row.names = rownames(tData[[1]]),
                                Cluster = UMAP_clustres,
                                stringsAsFactors = FALSE)
        Clusters = max(UMAP_clustres)
        CoordData = as.data.frame(Concat_UMAP,
                                  row.names = rownames(tData[[1]]))
        self$Multi_Omics$ClusterRes$GAUDI$HDBSCAN[[paste0("Clusters_", Clusters)]] = ClusterRes
        self$Multi_Omics$Fit$GAUDI$HDBSCAN[[paste0("Clusters_", Clusters)]]$Single_Omics_UMAPs = Single_Omics_UMAPs
        self$Multi_Omics$Fit$GAUDI$HDBSCAN[[paste0("Clusters_", Clusters)]]$Concatenated_UMAP = Concat_UMAP
        self$Multi_Omics$Fit$GAUDI$HDBSCAN[[paste0("Clusters_", Clusters)]]$Cluster_Method = ClusterMethod
        self$Multi_Omics$CoordData$GAUDI$HDBSCAN[[paste0("Factors_", ConcatUMAPDimensions)]] = CoordData
        self$Multi_Omics$Fit$GAUDI$HDBSCAN[[paste0("Clusters_", Clusters)]]$Param = CurrentParam

      }else if(ClusterMethod == "DBSCAN"){
        UMAP_clustres = dbscan::dbscan(Concat_UMAP, eps = DBSCANeps, minPts = DBSCANminPts)[[1]]

        #save results
        ClusterRes = data.frame(row.names = rownames(tData[[1]]),
                                Cluster = UMAP_clustres,
                                stringsAsFactors = FALSE)
        Clusters = max(UMAP_clustres)
        CoordData = as.data.frame(Concat_UMAP,
                                  row.names = rownames(tData[[1]]))
        self$Multi_Omics$ClusterRes$GAUDI$DBSCAN[[paste0("Clusters_", Clusters)]] = ClusterRes
        self$Multi_Omics$Fit$GAUDI$DBSCAN[[paste0("Clusters_", Clusters)]]$Single_Omics_UMAPs = Single_Omics_UMAPs
        self$Multi_Omics$Fit$GAUDI$DBSCAN[[paste0("Clusters_", Clusters)]]$Concatenated_UMAP = Concat_UMAP
        self$Multi_Omics$Fit$GAUDI$DBSCAN[[paste0("Clusters_", Clusters)]]$Cluster_Method = ClusterMethod
        self$Multi_Omics$CoordData$GAUDI$DBSCAN[[paste0("Factors_", ConcatUMAPDimensions)]] = CoordData
        self$Multi_Omics$Fit$GAUDI$DBSCAN[[paste0("Clusters_", Clusters)]]$Param = CurrentParam

      }else if(ClusterMethod == "Kmeans"){
        for(i in KmeansClusters){
          UMAP_clustres = stats::kmeans(Concat_UMAP, centers = i)[[1]]

          #save results
          ClusterRes = data.frame(row.names = rownames(tData[[1]]),
                                  Cluster = UMAP_clustres,
                                  stringsAsFactors = FALSE)
          Clusters = max(UMAP_clustres)
          CoordData = as.data.frame(Concat_UMAP,
                                    row.names = rownames(tData[[1]]))
          self$Multi_Omics$ClusterRes$GAUDI$Kmeans[[paste0("Clusters_", Clusters)]] = ClusterRes
          self$Multi_Omics$Fit$GAUDI$Kmeans[[paste0("Clusters_", Clusters)]]$Single_Omics_UMAPs = Single_Omics_UMAPs
          self$Multi_Omics$Fit$GAUDI$Kmeans[[paste0("Clusters_", Clusters)]]$Concatenated_UMAP = Concat_UMAP
          self$Multi_Omics$Fit$GAUDI$Kmeans[[paste0("Clusters_", Clusters)]]$Cluster_Method = ClusterMethod
          self$Multi_Omics$CoordData$GAUDI$Kmeans[[paste0("Factors_", ConcatUMAPDimensions)]] = CoordData
          self$Multi_Omics$Fit$GAUDI$Kmeans[[paste0("Clusters_", Clusters)]]$Param = CurrentParam
        }
      }
    },

    #' @description
        #' Hierarchical clustering for the feature weights calculated by the
        #' multi-omics integration methods MoCluster, MCIA, jNMF, iNMF and MOFA.
        #' Integrated omics feature datasets from $Multi_Omics$FeatureRes are
        #' scaled and concatenated, and distance matrices are calculated per
        #' multi-omics integration method. Hierarchical clustering is then
        #' performed to cluster features from different omics types.
    #' @param nFactors Number of factors (components) for which the feature
    #' hierarchical clustering is to be performed. The multi-omics integration
    #' methods need to be performed with this number of components (MoCluster
    #' and MCIA), clusters (jNMF and iNMF) or factors (MOFA).
    #' @param Distance Distance metric to be used for the calculation of the
    #' feature distance matrices. Must be one of: `"euclidean"` (default),
    #' `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"` or `"minkowski"`.
    #' If `"minkowski"` is selected, argument `"MinkowskiPower"` needs to be
    #' included.
    #' @param MinkowskiPower Power of the Minkowski distance. Default is NULL.
    #' @param Linkage Agglomeration method to be used for the hierarchical
    #' clustering. Must be one of: `"ward.D"`, `"ward.D2"` (default), `"single"`,
    #' `"complete"`, `"average"`, `"mcquitty"`, `"median"` or `"centroid"`.
    #' @returns A list of results. Feature distance matrix, stored by
    #' multi-omics integration method in
    #' $Multi_Omics$Feature_HClust$DistanceMatrix. Feature hierarchical
    #' clustering trees, stored by multi-omics integration method in
    #' $Multi_Omics$Feature_HClust$HClust. Feature dendrogram, stored by
    #' multi-omics integration method in $Plots$Dendrogram$Features$Multi_Omics.
    #' @references   * Meng C (2025). mogsa: Multiple omics data integrative
    #'   clustering and gene set analysis. doi:10.18129/B9.bioc.mogsa, R package
    #'   version 1.42.0
    #'   * Meng C, Kuster B, Culhane A, Gholami AM (2013). “A multivariate
    #'   approach to the integration of multi-omics datasets.” BMC Bioinformatics.
    #'   * Tsuyuzaki, K., and Nikaido, I. (2024). nnTensor: Non-Negative Tensor
    #'   Decomposition. R package version 1.3.0
    #'   * Chalise, P., Raghavan, R., and Fridley, B. (2025). IntNMF: Integrative
    #'   Clustering of Multiple Genomic Dataset. R package version 1.3.0
    #'   * Argelaguet R, Velten B, Arnol D, Dietrich S, Zenz T, Marioni JC,
    #'   Buettner F, Huber W, Stegle O (2018). “Multi‐Omics Factor Analysis—a
    #'   framework for unsupervised integration of multi‐omics data sets.”
    #'   Molecular Systems Biology, 14.
    run_Multi_Omics_Feature_HClust = function(nFactors = NULL,
                                              Distance = "euclidean",
                                              MinkowskiPower = NULL,
                                              Linkage = "ward.D2"){
      Data = self$Multi_Omics$FeatureRes
      ScaleData = Feature_Weight_Scaling(FeatureData = Data)
      MOFeatureHClust = Multi_Omics_Feature_Clustering(ScaledFeaturesData = ScaleData,
                                                       nFactors = nFactors,
                                                       Distance = Distance,
                                                       MinkowskiPower = MinkowskiPower,
                                                       Linkage = Linkage)
      self$Multi_Omics$Feature_HClust$DistanceMatrix = MOFeatureHClust$DistMat
      self$Multi_Omics$Feature_HClust$HClustTree = MOFeatureHClust$HClust
      self$Plots$Dendrogram$Features$Multi_Omics = MOFeatureHClust$Dendrogram
    },

    #----------------------------------------------------------------ENSEMBLE METHODS------------

    #Ensemble sample clustering by Multidimensional Scaling (MDS)
    #' @description
        #' Create an ensemble clustering result by combining the clustering
        #' results from MoCluster, MCIA, jNMF, iNMF, LRAcluster, COCA, MOFA,
        #' SNF and GAUDI. Borrows the idea of consensus cluster, by creating a
        #' sample similarity matrix counting the number of times samples cluster
        #' together. Hierarchical clustering is then performed on this matrix.
        #' Provides a more robust clustering result when performed with
        #' multiple cluster assignment, which can be achieved by running the
        #' multi-omics integration methods with a range of `Clusters`.
    #' @param Clusters Number of clusters to create from the hierarchical
    #' clustering result.
    #' @param SNFDistance Which distance used to calculate SNF results to be
    #' included in the ensemble clustering. This is needed to reduce the impact
    #' of SNF on the ensemble clustering assignment. Must be one of:
    #' `"euclidean squared"` (default), `"euclidean"`, `"manhattan"`,
    #' `"minkowksi 0.25"`, `"minkowski 0.5"`, `"minkowski 3"` or `"minkowski 4"`.
    #' @param Distance Distance Distance metric to be used for the calculation of the
    #' feature distance matrices. Must be one of: `"euclidean"` (default),
    #' `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"` or `"minkowski"`.
    #' If `"minkowski"` is selected, argument `"MinkowskiPower"` needs to be
    #' included.
    #' @param MinkowskiPower Power of the Minkowski distance. Default is NULL.
    #' @param Linkage Agglomeration method to be used for the hierarchical
    #' clustering. Must be one of: `"ward.D"`, `"ward.D2"`, `"single"`,
    #' `"complete"` (default), `"average"`, `"mcquitty"`, `"median"` or
    #' `"centroid"`.
    #' @returns A list of results. Sample distance matrix is stored in
    #' $Ensemble$Samples$MDS$DistanceMatrix by cluster number. Sample hierarchical
    #' clustering tree is stored in $Ensemble$Samples$MDS$HClust by cluster number.
    #' Ensemble sample cluster assignment is stored in
    #' $Ensemble$Samples$MDS$ClusterRes by cluster number. Sample dendrogram is
    #' stored in $Plots$Dendrogram$Ensemble$MDS by cluster number.
    run_Ensemble_Sample_MDS = function(Clusters = NULL,
                                       SNFDistance = "euclidean squared",
                                       Distance = "euclidean",
                                       MinkowskiPower = NULL,
                                       Linkage = "ward.D2"){
      Ens_Clust_Methods = names(self$Multi_Omics$ClusterRes)
      ClusterResults = list()

      for(method in Ens_Clust_Methods){
        if(method == "SNF"){
          Ens_Clust_SNF_Distance = names(self$Multi_Omics$ClusterRes$SNF)
          for(distance in Ens_Clust_SNF_Distance){
            if(distance == SNFDistance){
              ClusterResults$SNF[[distance]] = self$Multi_Omics$ClusterRes$SNF[[distance]]
            }
          }
        }else{
          ClusterResults[[method]] = self$Multi_Omics$ClusterRes[[method]]
        }
      }

      Ensemble_Data = Ensemble_Cluster_Data(Data = ClusterResults)
      Ensemble_Sample_IDs = rownames(Ensemble_Data[[1]])

      Ensemble_Clustering = Ensemble_Cluster(MethodClusterResults = Ensemble_Data,
                                             Sample_IDs = Ensemble_Sample_IDs,
                                             Distance = Distance,
                                             MinkowskiPower = MinkowskiPower,
                                             Linkage = Linkage,
                                             Clusters = Clusters)
      ClusterRes = data.frame(row.names = Ensemble_Sample_IDs,
                              Cluster = Ensemble_Clustering$ClusterRes,
                              stringsAsFactors = FALSE)
      EnsembleDendro = stats::as.dendrogram(Ensemble_Clustering$HClustRes)
      EnsembleClustnum = paste0("Clusters_", Clusters)
      self$Ensemble$Samples$MDS$ClusterRes[[EnsembleClustnum]] = ClusterRes
      self$Ensemble$Samples$MDS$DistanceMatrix[[EnsembleClustnum]] = Ensemble_Clustering$DistMat
      self$Ensemble$Samples$MDS$HClustTree[[EnsembleClustnum]] = Ensemble_Clustering$HClustRes
      self$Plots$Dendrogram$Ensemble$MDS[[EnsembleClustnum]] = EnsembleDendro
    },

    #Ensemble sample clustering by Canonical Correlation Analysis (CCA)
    #' @description
    #' Create an ensemble clustering result by aligning the embeddings results
    #' from MoCluster, MCIA, jNMF, iNMF, LRAcluster, COCA, MOFA, SNF and GAUDI.
    #' Uses CCA to align the embeddings from the multi-omics integration methods
    #' into a common latend space. For COCA and SNF, low-dimensional embeddings
    #' of samples are not provided, and one-hot encoding is used in place of
    #' embedding. One method is selected as reference embedding, with all other
    #' embeddings sequentially aligned to this reference using pairwise CCA. The
    #' aligned embedding is then used for k-means clustering.
    #' @param Clusters Number of clusters to create from the CCA aligned
    #' embedding.
    #' @param SNFDistance Which distance used to calculate SNF results to be
    #' included in the ensemble clustering. This is needed to reduce the impact
    #' of SNF on the ensemble clustering assignment. Must be one of:
    #' `"euclidean squared"` (default), `"euclidean"`, `"manhattan"`,
    #' `"minkowksi 0.25"`, `"minkowski 0.5"`, `"minkowski 3"` or `"minkowski 4"`.
    #' @param Reference Multi-omics integration embedding to be used as
    #' reference for CCA alignment.
    #' @param Factors Number of factors for which the multi-omics integration
    #' embedding is to be selected.
    #' @param GAUDIMethod Which clustering method is used to calculate GAUDI
    #' results to be included in the ensemble clustering. This is needed to
    #' reduce the impact of GAUDI on the ensemble clustering assignment. Must be
    #' one of: `"Kmeans"` (default), `"HDBSCAN"` or `"DBSCAN"`.
    #' @returns A list of results. CCA aligned embedding is stored in
    #' $Ensemble$Samples$CCA$Embedding by factor number. Ensemble sample cluster
    #' assignment is stored in $Ensemble$Samples$CCA$ClusterRes by factor and
    #' cluster number.
    run_Ensemble_Sample_CCA = function(Clusters = NULL,
                                       SNFDistance = "euclidean squared",
                                       Reference = "MoCluster",
                                       Factors = NULL,
                                       GAUDIMethod = "Kmeans"){
      #Data selection
      if(!is.null(self$Preprocessed_Omics$Matched_Data)){
        Data = self$Preprocessed_Omics$Matched_Data
      } else if(!is.null(self$Preprocessed_Omics$Data)){
        Data = self$Preprocessed_Omics$Data
      } else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      Samples = colnames(Data[[1]])

      #extract embeddings
      FactorNum = paste0("Factors_", Factors)
      EmbedList = list()
      #with CoordData
      Coord_Methods = names(self$Multi_Omics$CoordData)
      for(method in Coord_Methods){
        if(method == "GAUDI"){
          Embed = self$Multi_Omics$CoordData[[method]][[GAUDIMethod]][[FactorNum]]
        }else{
          Embed = self$Multi_Omics$CoordData[[method]][[FactorNum]]
        }
        if(!is.null(Embed)){
          Embed = Embed[match(Samples, rownames(Embed)), , drop = FALSE]
          colnames(Embed) <- paste0(method, "_Factor_", 1:ncol(Embed))
          EmbedList[[method]] = as.matrix(Embed)
        }else{
          warning(paste0("No embedding from ", method, " with ", Factors, " factors available, rerun ", method, " to include this method in CCA"))
        }
      }

      #with jNMF/iNMF
      ClustNum = paste0("Clusters_", Factors)
      NMFMethod = c("jNMF", "iNMF")
      for(method in NMFMethod){
        FitW = self$Multi_Omics$Fit[[method]][[ClustNum]]$W
        if(!is.null(dim(FitW))){
          FitW = FitW[match(Samples, rownames(FitW)), , drop = FALSE]
          colnames(FitW) = paste0(method, "_Factor_", 1:ncol(FitW))
          EmbedList[[method]] = as.matrix(FitW)
        }else{
          warning(paste0("No embedding from ", method, " with ", Factors, " factors available, rerun ", method, " to include this method in CCA"))
        }
      }

      #no coord data
      OneHotMethod = c("COCA", "SNF")
      for(method in OneHotMethod) {
        if(method == "SNF"){
          ClustRes = self$Multi_Omics$ClusterRes[[method]][[SNFDistance]][[ClustNum]]
        }else{
          ClustRes = self$Multi_Omics$ClusterRes[[method]][[ClustNum]]
        }
        if (!is.null(ClustRes)) {
          # Match sample order
          ClustRes = ClustRes[match(Samples, rownames(ClustRes)), , drop = FALSE]
          ClusterLabel = as.numeric(ClustRes$Cluster)

          # One-hot encode
          n_samples = length(ClusterLabel)
          n_clusters = Factors
          OneHot = matrix(0, nrow = n_samples, ncol = n_clusters)

          for (i in 1:n_samples) {
            if (!is.na(ClusterLabel[i]) && ClusterLabel[i] <= n_clusters) {
              OneHot[i, ClusterLabel[i]] = 1
            }
          }

          rownames(OneHot) = Samples
          colnames(OneHot) = paste0(method, "_Cluster_", 1:n_clusters)
          EmbedList[[method]] = as.matrix(OneHot)
        }else{
          warning(paste0("No embedding from ", method, " with ", Factors, " factors available, rerun ", method, " to include this method in CCA"))
        }
      }

      if(length(EmbedList) < 2) {
        stop(paste0("Fewer than 2 embeddings with ", Factors, " factors available, run additional methods to perform ensemble CCA"))
      }

      if(Reference %in% names(EmbedList)){
        ref_id = which(names(EmbedList) == Reference)
        EmbedList = c(EmbedList[ref_id], EmbedList[-ref_id])
      }else{
        stop(paste0("No ", Reference, " embeddings available with ", Factors, "factors, run ", Reference, " or change Reference variable"))
      }

      Ref = as.matrix(EmbedList[[1]])
      n_samples = nrow(Ref)

      Aligned = list()
      Aligned[[1]] = Ref

      for(i in 2:length(EmbedList)){
        CurrentEmbed = as.matrix(EmbedList[[i]])
        if(nrow(CurrentEmbed) != n_samples){
          Aligned[[i]] = CurrentEmbed
          next
        }

        tryCatch({
          n_CCA = min(ncol(Ref), ncol(CurrentEmbed), n_samples - 1)

          CCARes = stats::cancor(x = as.matrix(Ref[, 1:n_CCA, drop = FALSE]),
                                 y = as.matrix(CurrentEmbed[, 1:n_CCA, drop = FALSE]))

          if(!is.null(CCARes$ycoef) && is.matrix(CCARes$ycoef)){
            n_Canon = min(ncol(CCARes$ycoef), n_CCA)
            AlignedCurr = as.matrix(CurrentEmbed[, 1:n_CCA, drop = FALSE]) %*%
              as.matrix(CCARes$ycoef[, 1:n_Canon, drop = FALSE])
            Aligned[[i]] = AlignedCurr
          }else{
            Aligned[[i]] = CurrentEmbed
          }
        }, error = function(e) {
          Aligned[[i]] <<- CurrentEmbed
        })
      }

      Aligned_Embed = do.call(cbind, Aligned)
      rownames(Aligned_Embed) = Samples
      colnames(Aligned_Embed) = paste0("CCA_", 1:ncol(Aligned_Embed))

      Aligned_Clust = stats::kmeans(Aligned_Embed, centers = Clusters)[[1]]
      AlClusterRes = data.frame(row.names = Samples,
                                Cluster = Aligned_Clust,
                                stringsAsFactors = FALSE)

      EnsembleClustnum = paste0("Clusters_", Clusters)
      EnsembleFactnum = paste0("Embedding_", Factors)
      self$Ensemble$Samples$CCA[[EnsembleFactnum]]$ClusterRes[[EnsembleClustnum]] = AlClusterRes
      self$Ensemble$Samples$CCA[[EnsembleFactnum]]$EmbeddingMat = Aligned_Embed
    },

    #Feature hierarchical clustering
    #' @description
        #' Create an ensemble feature clustering result by combining the feature
        #' weight results from MoCluster, MCIA, jNMF, iNMF and MOFA. Feature
        #' weigths must be calculated for one number of components (MoCluster
        #' and MICA), clusters (jNMF and iNMF) and factors (MOFA) to be included.
        #' Two methods of combining the feature weights are available:
        #' `"Average"` and `"Concatenation"`. More details of the methods are
        #' available in the `param` `Method`.
    #' @param Method Method used for combining the feature weight matrices of
    #' the different multi-omics integration methods. Must be one of:
    #'    * `"Average"`: Feature weight matrices of the different omics types
    #'    are combined per multi-omics integration method. Feature distances are
    #'    calculated per method. The average feature distance between the five
    #'    methods is calculated. Hierachical clustering is performed on the
    #'    average feature distance matrix. Default method.
    #'    * `"Concatenation"`: Feature weight matrices of the different omics
    #'    types are combined per multi-omics intergration method, which are then
    #'    concatenated into a large matrix. A feature distance matrix as well as
    #'    a factor distance matrix is calculated. Hierarchical clustering is
    #'    then performed for both the features and the factors per method.
    #' @param nFactors Number of factors (components) for which the feature
    #' hierarchical clustering is to be performed. The multi-omics integration
    #' methods need to be performed with this number of components (MoCluster
    #' and MCIA), clusters (jNMF and iNMF) or factors (MOFA).
    #' @param Distance Distance Distance metric to be used for the calculation of the
    #' feature distance matrices. Must be one of: `"euclidean"` (default),
    #' `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"` or `"minkowski"`.
    #' If `"minkowski"` is selected, argument `"MinkowskiPower"` needs to be
    #' included.
    #' @param MinkowskiPower Power of the Minkowski distance. Default is NULL.
    #' @param Linkage Agglomeration method to be used for the hierarchical
    #' clustering. Must be one of: `"ward.D"`, `"ward.D2"`, `"single"`,
    #' `"complete"` (default), `"average"`, `"mcquitty"`, `"median"` or
    #' `"centroid"`.
    #' @returns A list of results. Feature distance matrix, stored by feature
    #' ensemble method in $Ensemble$Features[[`Method`]]$DistanceMatrix.
    #' Feature hierachical clustering tree, stored by feature ensemble method in
    #' $Ensemble$Features[[`Method`]]$HClust. Feature dendrogram, stored by
    #' feature ensemble method in $Plots$Dendrogram$Features. Factor distance
    #' matrix, stored in $Ensemble$Factors$Concatenation$DistanceMatrix. Only
    #' available for `Method` `"Concatenation"`. Factor hierarchical clustering
    #' tree, stored in $Ensemble$Factors$Concatenation$HClust. Only available
    #' for `Method` `"Concatenation"`. Factor dendrogram, stored in
    #' $Plots$Dendrogram$Factors$Concatenation. Only available for `Method`
    #' `"Concatenation"`.
    run_Ensemble_Integration_Feature = function(Method = "Average",
                                                nFactors = NULL,
                                                Distance = "euclidean",
                                                MinkowskiPower = NULL,
                                                Linkage = "ward.D2"
    ){
      Data = self$Multi_Omics$FeatureRes
      ScaleData = Feature_Weight_Scaling(FeatureData = Data)
      if(Method == "Average" | Method == "average" | Method == "Av" | Method == "av"){
        FeatureHClust = Feature_Method_Av_HClust(ScaledFeaturesData = ScaleData,
                                                 nFactors = nFactors,
                                                 Distance = Distance,
                                                 MinkowskiPower = MinkowskiPower,
                                                 Linkage = Linkage)
        self$Ensemble$Features$Average$DistanceMatrix = FeatureHClust$DistMat
        self$Ensemble$Features$Average$HClustTree = FeatureHClust$HClust
        self$Plots$Dendrogram$Features$Average = FeatureHClust$Dendrogram
      }else if(Method == "Concatenation" | Method == "concatenation" | Method == "Concat" | Method == "concat" | Method == "Con" | Method == "con"){
        FullFeatures = Concat_Feature_Matrix(ScaledFeaturesData = ScaleData,
                                             nFactors = nFactors)
        FeatureHClust = Feature_Concat_HClust(FullFeatureData = FullFeatures,
                                              Distance = Distance,
                                              MinkowskiPower = MinkowskiPower,
                                              Linkage = Linkage)
        self$Ensemble$Features$Concatenation$DistanceMatrix = FeatureHClust$FeatureDist
        self$Ensemble$Features$Concatenation$HClustTree = FeatureHClust$FeatureHClust
        self$Plots$Dendrogram$Features$Concatenation = FeatureHClust$FeatureDendro
        self$Ensemble$Factors$Concatenation$DistanceMatrix = FeatureHClust$FactorDist
        self$Ensemble$Factors$Concatenation$HClustTree = FeatureHClust$FactorHClust
        self$Plots$Dendrogram$Factors$Concatenation = FeatureHClust$FactorDendro
      }else{
        stop("Method needs to be either Average or Concatenation")
      }
    },

    #feature optimal cluster number
    #' @description
        #' Calculate the optimal number of feature clusters using the feature
        #' hierarchical clustering tree calculated in
        #' `run_Feature_Ensemble_HClust`. Will calculate the optimal number of
        #' feature clusters for both `"Average"` and `"Concatenation"` methods
        #' at the same time. Uses the as.clustrange function from the
        #' WeightedCluster package to calculate the optimal number of clusters.
    #' @param WeightedClusterStat The statistics included to calculate the
    #' optimal cluster number. Must be one of: `"all"` (default), `"noCH"` (all
    #' statistics except `"CH"` and `"CHsq"`), `"PBC"`, `"HG"`, `"HGSD"`,
    #' `"ASW"`, `"ASWw"`, `"CH"`, `"R2"`, `"CHsq"`, `"R2sq"` or `"HC"`.
    #' @param MaxClusters The maximum number of feature clusters to be
    #' considered by the statistics. Default is 20.
    #' @returns The optimal number of feature clusters, stored in
    #' $Ensemble$Features[[`Method`]]$OptimalClusters by
    #' `run_Feature_Ensemble_HClust` `Method` as available.
    #' @references Studer M (2013). “WeightedCluster Library Manual: A practical
    #' guide to creating typologies of trajectories in the social sciences with
    #' R.” LIVES Working Papers 24.
    run_Feature_Optimal_Cluster = function(WeightedClusterStat = "all",
                                           MaxClusters = 20){
      if(!is.null(self$Ensemble$Features$Average$HClust)){
        FeatureOptClust = Feature_Optimal_Clusters(HClust = self$Ensemble$Features$Average$HClust,
                                                   Dist = self$Ensemble$Features$Average$DistanceMatrix,
                                                   WeightedClusterStat = WeightedClusterStat,
                                                   MaxClusters = MaxClusters)
        self$Ensemble$Features$Average$OptimalClusters = FeatureOptClust$OptClust$K
      }
      if(!is.null(self$Ensemble$Features$Concatenation$HClust)){
        FeatureOptClust = Feature_Optimal_Clusters(HClust = self$Ensemble$Features$Concatenation$HClust,
                                                   Dist = self$Ensemble$Features$Concatenation$DistanceMatrix,
                                                   WeightedClusterStat = WeightedClusterStat,
                                                   MaxClusters = MaxClusters)
        self$Ensemble$Features$Concatenation$OptimalClusters = FeatureOptClust$OptClust$K
      }
      if(is.null(self$Ensemble$Features$Average$HClust) & is.null(self$Ensemble$Features$Concatenation$HClust)){
        stop("No feature distance matrix available, run run_Feature_Ensemble_HClust first")
      }
    },

    #factor optimal cluster number
    #' @description
        #' Calculate the optimal number of factor clusters using the factor
        #' hierarchical clustering tree calculated in
        #' `run_Feature_Ensemble_HClust` with `Method` `"Concatenation"`. Uses the
        #' as.clustrange function from the WeightedCluster package to calculate the
        #' optimal number of clusters.
    #' @param WeightedClusterStat The statistics included to calculate the
    #' optimal cluster number. Must be one of: `"all"` (default), `"noCH"` (all
    #' statistics except `"CH"` and `"CHsq"`), `"PBC"`, `"HG"`, `"HGSD"`,
    #' `"ASW"`, `"ASWw"`, `"CH"`, `"R2"`, `"CHsq"`, `"R2sq"` or `"HC"`.
    #' @param MaxClusters The maximum number of factor clusters to be
    #' considered by the statistics. Default is 10.
    #' @returns The optimal number of factor clusters, stored in
    #' $Ensemble$Factors$Concatenation$OptimalClusters.
    #' @references Studer M (2013). “WeightedCluster Library Manual: A practical
    #' guide to creating typologies of trajectories in the social sciences with
    #' R.” LIVES Working Papers 24.
    run_Factor_Optimal_Cluster = function(WeightedClusterStat = "all",
                                          MaxClusters = 10){
      if(!is.null(self$Ensemble$Factors$Concatenation$HClust)){
        FactorOptClust = Feature_Optimal_Clusters(HClust = self$Ensemble$Factors$Concatenation$HClust,
                                                  Dist = self$Ensemble$Factors$Concatenation$DistanceMatrix,
                                                  WeightedClusterStat = WeightedClusterStat,
                                                  MaxClusters = MaxClusters)
        self$Ensemble$Factors$OptimalClusters = FactorOptClust$OptClust$K
      }else{
        stop("No factor distance matrix available, run run_Feature_Ensemble_HClust with Method = Concatenation first")
      }
    },

    #Clustering of the feature HClust
    #' @description
        #' Cut the feature hierarchical clustering tree to create feature
        #' clusters. Can be done either using the calculated optimal feature
        #' cluster numbers, or with user argument.
    #' @param Clusters Number of clusters in which to cut the feature
    #' hierarchical clustering tree. Default is `NULL.` If `NULL`, the
    #' optimal feature cluster numbers calculated with
    #' `run_Feature_Optimal_Cluster` are used, and hierarchical clustering trees
    #' for all available `run_Feature_Ensemble_HClust` `Method` results. If an
    #' integer, the `run_Feature_Ensemble_HClust` `Method` whose hierarchical
    #' clustering tree is to be cut must be provided in parameter `Method`.
    #' @param Method Which hierarchical clustering tree must be cut. Only
    #' required if `Clusters` is an integer. Must be one of: `"Average"` or
    #' `"Concatenation"`.
    #' @returns A list of results. Feature clustering result, stored by
    #' `run_Feature_Ensemble_HClust` `Method` in
    #' $Ensemble$Features[[`Method`]]$ClusterRes. Feature dendrogram with
    #' cluster assignment, stored by `run_Feature_Ensemble_HClust` `Method` in
    #' $Plots$Dendrogram$Features.
    run_Feature_Clustering = function(Clusters = NULL, #Manual cluster number
                                      Method = NULL #Average or Concatenation
    ){
      if(is.null(Clusters)){
        if(!is.null(self$Ensemble$Features$Average$HClust)){
          Clusters = self$Ensemble$Features$Average$OptimalClusters
          FeatClust = Feature_Dendro_Clustering(FeatHClust = self$Ensemble$Features$Average$HClust,
                                                Clusters = Clusters)
          self$Ensemble$Features$Average$ClusterRes = FeatClust$Clusters
          self$Plots$Dendrogram$Features$Average = FeatClust$Dendro
        }
        if(!is.null(self$Ensemble$Features$Concatenation$HClust)){
          Clusters = self$Ensemble$Features$Concatenation$OptimalClusters
          FeatClust = Feature_Dendro_Clustering(FeatHClust = self$Ensemble$Features$Concatenation$HClust,
                                                Clusters = Clusters)
          self$Ensemble$Features$Concatenation$ClusterRes = FeatClust$Clusters
          self$Plots$Dendrogram$Features$Concatenation = FeatClust$Dendro
        }
        if(is.null(self$Ensemble$Features$Average$HClust) & is.null(self$Ensemble$Features$Concatenation$HClust)){
          stop("No feature hclust available, run run_Feature_Ensemble_HClust first")
        }
      }else{
        if(Method == "Average" | Method == "average" | Method == "Av" | Method == "av"){
          if(!is.null(self$Ensemble$Features$Average$HClust)){
            FeatClust = Feature_Dendro_Clustering(FeatHClust = self$Ensemble$Features$Average$HClust,
                                                  Clusters = Clusters)
            self$Ensemble$Features$Average$ClusterRes = FeatClust$Clusters
            self$Plots$Dendrogram$Features$Average = FeatClust$Dendro
          }else{
            stop("No feature hclust available, run run_Feature_Ensemble_HClust first")
          }
        }else if(Method == "Concatenation" | Method == "concatenation" | Method == "Concat" | Method == "concat" | Method == "Con" | Method == "con"){
          if(!is.null(self$Ensemble$Features$Concatenation$HClust)){
            FeatClust = Feature_Dendro_Clustering(FeatHClust = self$Ensemble$Features$Concatenation$HClust,
                                                  Clusters = Clusters)
            self$Ensemble$Features$Concatenation$ClusterRes = FeatClust$Clusters
            self$Plots$Dendrogram$Features$Concatenation = FeatClust$Dendro
          }else{
            stop("No feature hclust available, run run_Feature_Ensemble_HClust first")
          }
        }else{
          stop("Unknown feature ensemble method provided, please select either Average or Concatenation")
        }
      }
    },

    #Clustering of the factor HClust
    #' @description
        #' Cut the factor hierarchical clustering tree to create factor
        #' clusters. Can be done either using the calculated optimal factor
        #' cluster numbers, or with user argument.
    #' @param Clusters Number of clusters in which to cut the factor
    #' hierarchical clustering tree. Default is `NULL.` If `NULL`, the
    #' optimal factor cluster numbers calculated with
    #' `run_Factor_Optimal_Cluster` are used.
    #' @returns A list of results. Factor clustering result, stored in
    #' $Ensemble$Features$Concatenation$ClusterRes. Factor dendrogram with
    #' cluster assignment, stored in $Plots$Dendrogram$Factors$Concatenation.
    run_Factor_Clustering = function(Clusters = NULL){
      if(!is.null(self$Ensemble$Factors$Concatenation$HClust)){
        if(is.null(Clusters)){
          Clusters = self$Ensemble$Factors$Concatenation$OptimalClusters
        }
        FactClust = Feature_Dendro_Clustering(FeatHClust = self$Ensemble$Factors$Concatenation$HClust,
                                              Clusters = Clusters)
        self$Ensemble$Factors$Concatenation$ClusterRes = FactClust$Clusters
        self$Plots$Dendrogram$Factors$Concatenation = FactClust$Dendro
      }
    },

    #overrepresentation analysis
    #' @description
        #' Perform overrepresentation analysis from package ClusterProfiler
        #' on the clustered ensembled features to determine if features with
        #' similar biological functions group together in across multi-omics
        #' integration methods. Uses enrichment GO categories from the genome
        #' wide annotation database for human as provided in package
        #' org.Hs.eg.db. Requires ensemble feature hierarchical clustering to
        #' have been performed with `run_Feature_Ensemble_HClust` and a cluster
        #' assignment to have been created with `run_Feature_Clustering`. Only
        #' available for gene-based omics types (genomics, transcriptomics,
        #' methylation etc.).
    #' @param Method Which ensemble feature clustering result is to be used.
    #' Must be one of" `"Average"` or `"Concatenation"`. Default is `"Average"`.
    #' @param OmicsName For which omics dataset the overrepresentation analysis
    #' is to be performed. Only genetics-based omics types can be utilized in the
    #' genome wide annotation database.
    #' @param GeneNameType Keytype of the feature name for the omics dataset.
    #' Must be one of: `"ENTREZID"`, `"PFAM"`, `"IPI"`, `"PROSITE"`, `"ACCNUM"`,
    #' `"ALIAS"`, `"CHR"`, `"CHRLOC"`, `"CHRLOCEND"`, `"ENZYME"`, `"MAP"`,
    #' `"PATH"`, `"PMID"`, `"REFSEQ"`, `"SYMBOL"`, `"UNIGENE"`, `"ENSEMBL"`,
    #' `"ENSEMBLPROT"`, `"ENSEMBLTRANS"`, `"GENENAME"`, `"UNIPROT"`, `"GO"`,
    #' `"EVIDENCE"`, `"ONTOLOGY"`, `"GOALL"`, `"EVIDENCEALL"`, `"ONTOLOGYALL"`,
    #' `"OMIM"` or `"UCSCKG"`.
    #' @param Subontologies Which ontology to use. Must be one of: `"ALL"`
    #' (default), `"MF"`, `"BP"` or `"CC"`.
    #' @param pValue Cutoff value for the p-value. Default is 0.05.
    #' @param pAdjustment Method for the p-value adjustment. Must be one of:
    #' `"BH"` (default), `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`,
    #' `"BY"`, `"fdr"` or `"none"`.
    #' @param qValue Cutoff for the q-value. Default is 0.2.
    #' @param MinGenesPerTerm Minimal number of genes annotated per ontology
    #' term to be included for testing. Default is 10.
    #' @param MaxGenesPerTerm Maximum number of genes annotated per ontology
    #' term to be included for testing. Default is 500.
    #' @returns A list of results. A dataframe with all features of the omics
    #' type per cluster of ensembled features. Stored by `Method` in
    #' $Ensemble$Features[[`Method`]]$ORA[[`OmicsName`]][[`Clusters`]]$FeatureList.
    #' Enrichment result of class `enrichResult` per cluster of ensembled
    #' features, if that cluster contains features of the selected omics type
    #' and enrichment analysis resulted in enriched terms. Stored by `Method` in
    #' $Ensemble$Features[[`Method`]]$ORA[[`OmicsName`]][[`Clusters`]]$Enrichment.
    #' @references   * Yu G (2024). “Thirteen years of clusterProfiler.” The
    #'   Innovation, 5(6), 100722.
    #'   * Xu S, Hu E, Cai Y, Xie Z, Luo X, Zhan L, Tang W, Wang Q, Liu B, Wang R,
    #'   Xie W, Wu T, Xie L, Yu G (2024). “Using clusterProfiler to characterize
    #'   multiomics data.” Nature Protocols, 19(11).
    #'   * Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, Feng T, Zhou L, Tang W, Zhan L,
    #'   Fu x, Liu S, Bo X, Yu G (2021). “clusterProfiler 4.0: A universal
    #'   enrichment tool for interpreting omics data.” The Innovation, 2(3).
    #'   * Yu G, Wang L, Han Y, He Q (2012). “clusterProfiler: an R package for
    #'   comparing biological themes among gene clusters.” OMICS: A Journal of
    #'   Integrative Biology, 16(5), 284-287.
    #'   * Carlson M (2025). "Genome wide annotation for Human, primarily based on
    #'   mapping using Entrez Gene identifiers." Bioconductor.
    run_Over_Represenation = function(Method = "Average", #Feature hierarchical clustering method, either Average or Concatenation
                                      OmicsName,
                                      GeneNameType,
                                      Subontologies = "ALL",
                                      pValue = 0.05,
                                      pAdjustment = "BH",
                                      qValue = 0.2,
                                      MinGenesPerTerm = 10,
                                      MaxGenesPerTerm = 500
    ){
      if(Method == "Average" | Method == "average" | Method == "Av" | Method == "av"){
        FeatureData = self$Ensemble$Features$Average$ClusterRes
      }else if(Method == "Concatenation" | Method == "concatenation" | Method == "Concat" | Method == "concat" | Method == "Con" | Method == "con"){
        FeatureData = self$Ensemble$Features$Concatenation$ClusterRes
      }else{
        stop("Unknown method of feature hierarchical clustering, please select either Average or Concatenation")
      }

      EnrichResults = list()
      FeatureClusData = tibble::rownames_to_column(FeatureData, var = "FeatureOmics")
      FeatureClusData = transform(FeatureClusData, Feature = sub("(.*)_.*", "\\1", FeatureOmics), Omics = sub(".*_", "", FeatureOmics))
      OmicsFeatureClusData = subset(FeatureClusData, Omics == OmicsName)
      FullFeatureList = as.vector(OmicsFeatureClusData$Feature)
      MaxClust = max(OmicsFeatureClusData$Cluster)
      for(i in 1:MaxClust){
        OmicsFeatureClusi = subset(OmicsFeatureClusData, Cluster == i)
        ClusterFeatureList = as.vector(OmicsFeatureClusi$Feature)
        ClusterName = paste0("Cluster_", i)
        ClusterEnrich = ORA_GO(FullList = FullFeatureList,
                               SelectList = ClusterFeatureList,
                               GeneNameType = GeneNameType,
                               Subontologies = Subontologies,
                               pValue = pValue,
                               pAdjustment = pAdjustment,
                               qValue = qValue,
                               MinGenesPerTerm = MinGenesPerTerm,
                               MaxGenesPerTerm = MaxGenesPerTerm)
        EnrichResults[[ClusterName]]$FeatureList = ClusterFeatureList
        EnrichResults[[ClusterName]]$Enrichment = ClusterEnrich
      }

      if(Method == "Average" | Method == "average" | Method == "Av" | Method == "av"){
        self$Ensemble$Features$Average$ORA[[OmicsName]] = EnrichResults
      }else if(Method == "Concatenation" | Method == "concatenation" | Method == "Concat" | Method == "concat" | Method == "Con" | Method == "con"){
        self$Ensemble$Features$Concatenation$ORA[[OmicsName]] = EnrichResults
      }
    },

    #----------------------------------------------------------------PLOTTING METHODS-------------

    #heatmap
    #' @description
        #' Interactive heatmap with sample clustering results per omics type.
    #' @param FeatureClustering Whether to order the features on the Y-axis
    #' based on hierarchical clustering. Default is `TRUE`. If `FALSE`, the
    #' features will be ordered as in the preprocessed datasets.
    #' @param MethodResults For which clustering assignment the heatmap is to be
    #' plotted. Must be one of: `"all"` (default), `"MoCluster"`, `"jNMF"`, `"iNMF"`,
    #' `"LRAcluster"`, `"COCA"`, `"GAUDI"`, `"SNF"`, `"MDS"` or `"CCA"`. If all is
    #' selected, heatmaps with cluster assignments are created for all run
    #' methods.
    #' @param YFontSize Font size of the feature labels on the Y-axis. Default
    #' is 5.
    #' @param XFontSize Font size of the sample labels on the X-axis. Default
    #' is 5.
    #' @returns A list of plots. Heatmap with cluster assignment for each omics
    #' type per multi-omics or ensemble integration clustering results. Stored
    #' by integration method and cluster number for each omics type in
    #' $Plots$Cluster_Heatmap.
    plot_Heatmap_Cluster = function(FeatureClustering = TRUE, #whether to order the features based on hclust
                                    MethodResults = "all", #which multi-omics method to be plotted (method name, ensemble or all)
                                    YFontSize = 5,
                                    XFontSize = 5
    ){
      if(!is.null(self$Preprocessed_Omics$Matched_Data)){
        Data = self$Preprocessed_Omics$Matched_Data
      }else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      #per multi-omics integration method heatmap
      HeatmapClusterMethods = names(self$Multi_Omics$ClusterRes)
      HeatmapOmics = names(Data)

      for(method in HeatmapClusterMethods){
        if(method == "SNF"){
          if(MethodResults == "SNF" | MethodResults == "all"){
            HeatmapSNFDistances = names(self$Multi_Omics$ClusterRes$SNF)
            for(distance in HeatmapSNFDistances){
              HeatmapSNFClustNum = names(self$Multi_Omics$ClusterRes$SNF[[distance]])
              for(clustnum in HeatmapSNFClustNum){
                ClustData = self$Multi_Omics$ClusterRes$SNF[[distance]][[clustnum]]
                n = strsplit(clustnum, split = "_")[[1]][2]
                for(omics in HeatmapOmics){
                  Title = paste0("Heatmap of ", omics, ", ", n, " clusters by ", method, " ", distance)
                  Plot = Heatmap_Cluster(Data = Data[[omics]],
                                         SampleClusters = ClustData,
                                         FeatureClusters = FeatureClustering,
                                         YFontSize = YFontSize,
                                         XFontSize = XFontSize,
                                         Title = Title)
                  self$Plots$Cluster_Heatmap$SNF[[distance]][[clustnum]][[omics]] = Plot
                }
              }
            }
          }else{
            next
          }
        }else if(method == "GAUDI"){
          if(MethodResults == "GAUDI" | MethodResults == "all"){
            HeatmapGAUDIClustMethod = names(self$Multi_Omics$ClusterRes$GAUDI)
            for(clustmethod in HeatmapGAUDIClustMethod){
              HeatmapGAUDIClustNum = names(self$Multi_Omics$ClusterRes$GAUDI[[clustmethod]])
              for(clustnum in HeatmapGAUDIClustNum){
                ClustData = self$Multi_Omics$ClusterRes$GAUDI[[clustmethod]][[clustnum]]
                n = strsplit(clustnum, split = "_")[[1]][2]
                for(omics in HeatmapOmics){
                  Title = paste0("Heatmap of ", omics, ", ", n, " clusters by ", method, " , clustered with ", clustmethod)
                  Plot = Heatmap_Cluster(Data = Data[[omics]],
                                         SampleClusters = ClustData,
                                         FeatureClusters = FeatureClustering,
                                         YFontSize = YFontSize,
                                         XFontSize = XFontSize,
                                         Title = Title)
                  self$Plots$Cluster_Heatmap$GAUDI[[clustmethod]][[clustnum]][[omics]] = Plot
                }
              }
            }
          }else{
            next
          }
        } else{
          if(MethodResults == method | MethodResults == "all"){
            HeatmapClustNum = names(self$Multi_Omics$ClusterRes[[method]])
            for(clustnum in HeatmapClustNum){
              ClustData = self$Multi_Omics$ClusterRes[[method]][[clustnum]]
              n = strsplit(clustnum, split = "_")[[1]][2]
              for(omics in HeatmapOmics){
                Title = paste0("Heatmap of ", omics, ", ", n, " clusters by ", method)
                Plot = Heatmap_Cluster(Data = Data[[omics]],
                                       SampleClusters = ClustData,
                                       FeatureClusters = FeatureClustering,
                                       YFontSize = YFontSize,
                                       XFontSize = XFontSize,
                                       Title = Title)
                self$Plots$Cluster_Heatmap[[method]][[clustnum]][[omics]] = Plot
              }
            }
          }else{
            next
          }
        }
      }

      #Heatmap for consensus clustering result
      if(!is.null(self$Ensemble$Samples$MDS$ClusterRes) & (MethodResults == "MDS" | MethodResults == "all")){
        ConsensusClustNum = names(self$Ensemble$Samples$MDS$ClusterRes)
        for(clustnum in ConsensusClustNum){
          ClustData = self$Ensemble$Samples$MDS$ClusterRes[[clustnum]]
          n = strsplit(clustnum, split = "_")[[1]][2]
          for(omics in HeatmapOmics){
            Title = paste0("Heatmap of ", omics, ", ", n, " clusters by MDS ensemble")
            Plot = Heatmap_Cluster(Data = Data[[omics]],
                                   SampleClusters = ClustData,
                                   FeatureClusters = FeatureClustering,
                                   YFontSize = YFontSize,
                                   XFontSize = XFontSize,
                                   Title = Title)
            self$Plots$Cluster_Heatmap$Ensemble$MDS[[clustnum]][[omics]] = Plot
          }
        }
      }
      if(!is.null(self$Ensemble$Samples$CCA) & (MethodResults == "CCA" | MethodResults == "all")){
        CCAFactNum = names(self$Ensemble$Samples$CCA)
        for(factnum in CCAFactNum){
          CCAClustNum = names(self$Ensemble$Samples$CCA[[factnum]]$ClusterRes)
          for(clustnum in CCAClustNum){
            ClustData = self$Ensemble$Samples$CCA[[factnum]]$ClusterRes[[clustnum]]
            n = strsplit(clustnum, split = "_")[[1]][2]
            for(omics in HeatmapOmics){
              Title = paste0("Heatmap of ", omics, ", ", n, " clusters by CCA ensemble")
              Plot = Heatmap_Cluster(Data = Data[[omics]],
                                     SampleClusters = ClustData,
                                     FeatureClusters = FeatureClustering,
                                     YFontSize = YFontSize,
                                     XFontSize = XFontSize,
                                     Title = Title)
              self$Plots$Cluster_Heatmap$Ensemble$CCA[[clustnum]][[omics]] = Plot
            }
          }
        }
      }

    },

    #multi omics heatmap
    #' @description
        #' Non-interactive heatmap with all omics heatmaps in a single plot as
        #' well as clustering results and metadata features.
    #' @param LegendNames Names for the legends of the different omics heatmaps.
    #' If `NULL`, the omics names as used in the EINS object are used. Otherwise,
    #' a string vector with the same length as number of omics datasets should be
    #' provided.
    #' @param MetadataFeatures Names of the sample metadata features to be
    #' displayed on the heatmaps. These metadata features need to be available
    #' in the metadata files for all omics types. A string vector with the names
    #' of the metadata features as in the metadata files should be provided.
    #' @returns A list of plots. Multiple omics heatmaps with cluster assignment
    #' and metadata features per multi-omics or ensemble integration clustering
    #' result. Stored by integration method and cluster number in
    #' $Plots$Multi_Omics_Heatmap
    plot_Multi_Omics_Heatmap = function(LegendNames = NULL,
                                        MetadataFeatures = NULL){
      #select correct data for heatmap
      if(!is.null(self$Preprocessed_Omics$Matched_Data)){
        Data = self$Preprocessed_Omics$Matched_Data
      }else if(!is.null(self$Omics$Sample_Data)){
        Data = self$Omics$Sample_Data
      }else{
        Data = self$Omics$Raw_Data
      }

      #select the metadata
      Metadata = self$Omics$Metadata

      #create heatmap without cluster annotation
      NoClusterHeatmap = MultiOmicsHeatmap(OmicsData = Data,
                                           LegendNames = LegendNames,
                                           MetadataColumn = MetadataFeatures,
                                           MetaData = Metadata)
      self$Plots$Multi_Omics_Heatmap$No_Cluster = NoClusterHeatmap

      #create heatmaps with cluster annotation for all methods
      HeatmapMethods = names(self$Multi_Omics$ClusterRes)

      for(method in HeatmapMethods){
        if(method == "SNF"){
          HeatmapDistances = names(self$Multi_Omics$ClusterRes$SNF)
          for(distance in HeatmapDistances){
            HeatmapClustNum = names(self$Multi_Omics$ClusterRes$SNF[[distance]])
            for(clust in HeatmapClustNum){
              SNFClustRes = self$Multi_Omics$ClusterRes$SNF[[distance]][[clust]]
              ClusterHeatmap = MultiOmicsHeatmapClustered(OmicsData = Data,
                                                          ClusterRes = SNFClustRes,
                                                          LegendNames = LegendNames,
                                                          MetadataColumn = MetadataFeatures,
                                                          MetaData = Metadata)
              self$Plots$Multi_Omics_Heatmap$SNF[[distance]][[clust]] = ClusterHeatmap$Heatmap
            }
          }
        }else if(method == "GAUDI"){
          HeatmapClustMethod = names(self$Multi_Omics$ClusterRes$GAUDI)
          for(clustmethod in HeatmapClustMethod){
            HeatmapClustNum = names(self$Multi_Omics$ClusterRes$GAUDI[[clustmethod]])
            for(clust in HeatmapClustNum){
              GAUDIClustRes = self$Multi_Omics$ClusterRes$GAUDI[[clustmethod]][[clust]]
              ClusterHeatmap = MultiOmicsHeatmapClustered(OmicsData = Data,
                                                          ClusterRes = GAUDIClustRes,
                                                          LegendNames = LegendNames,
                                                          MetadataColumn = MetadataFeatures,
                                                          MetaData = Metadata)
              self$Plots$Multi_Omics_Heatmap$GAUDI[[clustmethod]][[clust]] = ClusterHeatmap$Heatmap
            }
          }
        }else{
          HeatmapClustNum = names(self$Multi_Omics$ClusterRes[[method]])
          for(clust in HeatmapClustNum){
            ClustRes = self$Multi_Omics$ClusterRes[[method]][[clust]]
            ClusterHeatmap = MultiOmicsHeatmapClustered(OmicsData = Data,
                                                        ClusterRes = ClustRes,
                                                        LegendNames = LegendNames,
                                                        MetadataColumn = MetadataFeatures,
                                                        MetaData = Metadata)
            self$Plots$Multi_Omics_Heatmap[[method]][[clust]] = ClusterHeatmap$Heatmap
          }
        }
      }

      #ensemble clustering multi-omics heatmap
      if(!is.null(self$Ensemble$Samples$MDS$ClusterRes)){
        EnsembleClustNum = names(self$Ensemble$Samples$MDS$ClusterRes)
        for(clust in EnsembleClustNum){
          ClustRes = self$Ensemble$Samples$MDS$ClusterRes[[clust]]
          ClusterHeatmap = MultiOmicsHeatmapClustered(OmicsData = Data,
                                                      ClusterRes = ClustRes,
                                                      LegendNames = LegendNames,
                                                      MetadataColumn = MetadataFeatures,
                                                      MetaData = Metadata)
          self$Plots$Multi_Omics_Heatmap$Ensemble$MDS[[clust]] = ClusterHeatmap
        }
      }
      if(!is.null(self$Ensemble$Samples$CCA)){
        CCAFactNum = names(self$Ensemble$Samples$CCA)
        for(factnum in CCAFactNum){
          EnsembleClustNum = names(self$Ensemble$Samples$CCA[[factnum]]$ClusterRes)
          for(clust in EnsembleClustNum){
            ClustRes = self$Ensemble$Samples$CCA[[factnum]]$ClusterRes[[clust]]
            ClusterHeatmap = MultiOmicsHeatmapClustered(OmicsData = Data,
                                                        ClusterRes = ClustRes,
                                                        LegendNames = LegendNames,
                                                        MetadataColumn = MetadataFeatures,
                                                        MetaData = Metadata)
            self$Plots$Multi_Omics_Heatmap$Ensemble$CCA[[clust]] = ClusterHeatmap
          }
        }
      }
    },

    #scatterplot
    #' @description
        #' Interactive scatterplot of samples in the factors
        #' (components/dimensions) for methods which have CoordData (MoCluster,
        #' MCIA, LRAcluster, MOFA and GAUDI).
    #' @param XAxis Integer, factor (component/dimension) to plot on the X-axis.
    #' Default is 1.
    #' @param YAxis Integer, factor (component/dimension) to plot on the Y-axis.
    #' Default is 2.
    #' @param YFontSize Size of the Y-axis title. Default is 15.
    #' @param XFontSize Size of the X-axis title. Default is 15.
    #' @param Width Width of the scatterplot. Default is 5.
    #' @param Height Height of the scatterplot. Default is 5.
    #' @param MetadataFeature Name of the metadata feature used as the shape of
    #' the dots. Must be a string with the name of the metadata feature as in
    #' the metadata files. This metadata feature needs to be available in the
    #' metadatafiles of all omics types.
    #' @returns A list of plots. Interactive scatterplots, with dots colored by
    #' cluster assignment and shaped by `MetadataFeature`. Stored by integration
    #' method and cluster number in $Plots$Scatterplot
    plot_Scatterplot = function(XAxis = 1, #factor to plot on x-axis
                                YAxis = 2, #factor to plot on y-axis
                                YFontSize = 15, #size of y label
                                XFontSize = 15, #size of x label
                                Width = 5, #width of plot
                                Height = 5, #height of plot
                                MetadataFeature = NULL
    ){
      MetaData = self$Omics$Metadata
      #create vector with metadata
      MatchedMetaData = Match_Metadata(MetaData)
      ReqMetaData = MatchedMetaData[MetadataFeature]

      ScatterplotMethods = names(self$Multi_Omics$CoordData)
      for(method in ScatterplotMethods){
        if(method == "GAUDI"){
          ScatterplotGAUDIClustMethod = names(self$Multi_Omics$CoordData$GAUDI)
          for(clustmethod in ScatterplotGAUDIClustMethod){
            ScatterplotClustNum = names(self$Multi_Omics$CoordData$GAUDI[[clustmethod]])
            for(clustnum in ScatterplotClustNum){
              CoordData = cbind(self$Multi_Omics$CoordData$GAUDI[[clustmethod]][[clustnum]], Clusters = as.factor(self$Multi_Omics$ClusterRes$GAUDI[[clustmethod]][[clustnum]]$Cluster), Metadata = as.factor(ReqMetaData[[1]]))
              n = strsplit(clustnum, split = "_")[[1]][2]
              Title = paste0("Scatterplot from ", method, " with ", n, " clusters, clustered by ", clustmethod)
              Plot = Scatterplot(Data = CoordData,
                                 XAxis = XAxis,
                                 YAxis = YAxis,
                                 Width = Width,
                                 Height = Height,
                                 Title = Title)
              self$Plots$Scatterplot$GAUDI[[clustmethod]][[clustnum]] = Plot
            }
          }
        }else{
          ScatterplotClustNum = names(self$Multi_Omics$CoordData[[method]])
          for(clustnum in ScatterplotClustNum){
            if(!is.null(self$Multi_Omics$ClusterRes[[method]][[clustnum]])){
              CoordData = cbind(self$Multi_Omics$CoordData[[method]][[clustnum]], Clusters = as.factor(self$Multi_Omics$ClusterRes[[method]][[clustnum]]$Cluster), Metadata = as.factor(ReqMetaData[[1]]))
              n = strsplit(clustnum, split = "_")[[1]][2]
              Title = paste0("Scatterplot from ", method, " with ", n, " clusters")
            }else{
              CoordData = cbind(self$Multi_Omics$CoordData[[method]][[clustnum]], Metadata = as.factor(ReqMetaData[[1]]))
              Title = paste0("Scatterplot from ", method)
            }
            Plot = Scatterplot(Data = CoordData,
                               XAxis = XAxis,
                               YAxis = YAxis,
                               Width = Width,
                               Height = Height,
                               Title = Title)
            self$Plots$Scatterplot[[method]][[clustnum]] = Plot
          }
        }
      }
    },


    #Create cluster + metadata colored plot
    #' @description
        #' Table with sample cluster assignment and metadata features, with each
        #' column colored by level.
    #' @param MetadataFeatures Names of the sample metadata features to be
    #' included in the table. These metadata features need to be available
    #' in the metadata files for all omics types. A string vector with the names
    #' of the metadata features as in the metadata files should be provided.
    #' @returns A list of plots. Cluster assignment and metadata feature table,
    #' stored by integration method and cluster number in
    #' $Plots$Clusters_Metadata.
    plot_Clusters_Metadata = function(MetadataFeatures = NULL #Columns of metadata to keep
    ){
      ClustersMetadataMethods = names(self$Multi_Omics$ClusterRes)
      MetaData = Match_Metadata(MetadataList = self$Omics$Metadata)

      for(method in ClustersMetadataMethods){
        if(method == "SNF"){
          ClustersMetadataDistances = names(self$Multi_Omics$ClusterRes$SNF)
          for(distance in ClustersMetadataDistances){
            ClustersMetadataClustNum = names(self$Multi_Omics$ClusterRes$SNF[[distance]])
            for(clustnum in ClustersMetadataClustNum){
              ClustData = self$Multi_Omics$ClusterRes$SNF[[distance]][[clustnum]]
              Plot = Plot_Clusters_Metadata(MetadataColumns = MetadataFeatures,
                                            MatchedMetadata = MetaData,
                                            ClusterResults = ClustData)
              self$Plots$Clusters_Metadata$SNF[[distance]][[clustnum]] = Plot
            }
          }
        } else if(method == "GAUDI"){
          ClustersMetadataClustMethod = names(self$Multi_Omics$ClusterRes$GAUDI)
          for(clustmethod in ClustersMetadataClustMethod){
            ClustersMetadataClustNum = names(self$Multi_Omics$ClusterRes$GAUDI[[clustmethod]])
            for(clustnum in ClustersMetadataClustNum){
              ClustData = self$Multi_Omics$ClusterRes$GAUDI[[clustmethod]][[clustnum]]
              Plot = Plot_Clusters_Metadata(MetadataColumns = MetadataFeatures,
                                            MatchedMetadata = MetaData,
                                            ClusterResults = ClustData)
              self$Plots$Clusters_Metadata$GAUDI[[clustmethod]][[clustnum]] = Plot
            }
          }
        } else{
          ClustersMetadataNum = names(self$Multi_Omics$ClusterRes[[method]])
          for(clustnum in ClustersMetadataNum){
            ClustData = self$Multi_Omics$ClusterRes[[method]][[clustnum]]
            Plot = Plot_Clusters_Metadata(MetadataColumns = MetadataFeatures,
                                          MatchedMetadata = MetaData,
                                          ClusterResults = ClustData)
            self$Plots$Clusters_Metadata[[method]][[clustnum]] = Plot
          }
        }
      }

      #single omics clustering results
      if(!is.null(self$Single_Omics$ClusterRes)){
        ClustersMetadataSingleOmics = names(self$Single_Omics$ClusterRes)
        for(omics in ClustersMetadataSingleOmics){
          ClustData = self$Single_Omics$ClusterRes[[omics]]
          Plot = Plot_Clusters_Metadata(MetadataColumns = MetadataFeatures,
                                        MatchedMetadata = MetaData,
                                        ClusterResults = ClustData)
          self$Plots$Clusters_Metadata$Single_Omics[[omics]] = Plot
        }
      }

      #ensemble clustering results
      if(!is.null(self$Ensemble$Samples$MDS$ClusterRes)){
        ClustersMetadataEnsemble = names(self$Ensemble$Samples$MDS$ClusterRes)
        for(clustnum in ClustersMetadataEnsemble){
          ClustData = self$Ensemble$Samples$MDS$ClusterRes[[clustnum]]
          Plot = Plot_Clusters_Metadata(MetadataColumns = MetadataFeatures,
                                        MatchedMetadata = MetaData,
                                        ClusterResults = ClustData)
          self$Plots$Clusters_Metadata$Ensemble$MDS[[clustnum]] = Plot
        }
      }
      if(!is.null(self$Ensemble$Samples$CCA)){
        CCAFactNum = names(self$Ensemble$Samples$CCA)
        for(factnum in CCAFactNum){
          ClustersMetadataEnsemble = names(self$Ensemble$Samples$CCA[[factnum]]$ClusterRes)
          for(clustnum in ClustersMetadataEnsemble){
            ClustData = self$Ensemble$Samples$CCA[[factnum]]$ClusterRes[[clustnum]]
            Plot = Plot_Clusters_Metadata(MetadataColumns = MetadataFeatures,
                                          MatchedMetadata = MetaData,
                                          ClusterResults = ClustData)
            self$Plots$Clusters_Metadata$Ensemble$CCA[[clustnum]] = Plot
          }
        }
      }
    },

    #Create Sankey plot comparing number of clusters in same method
    #' @description
        #' Interactive Sankey plot of the sample cluster results per method.
        #' Compares the cluster assignments of a single method with multiple
        #' cluster numbers. Can only be created if multiple cluster numbers are
        #' calculated per multi-omics integration method.
    #' @param MetadataFeature String of the sample metadata feature for which the
    #' samples should be colored. The metadata feature needs to be available
    #' in the metadata files for all omics types.
    #' @returns A list of plots. Sankey plot with the results of multiple
    #' cluster assignments from a single multi-omics integration method. Stored
    #' by multi-omics integration method in $Plots$Sankey_Clusters.
    plot_Sankey_Clusters = function(MetadataFeature = NULL #Column of metadata used for coloring samples
    ){
      Methods = names(self$Multi_Omics$ClusterRes)
      MetaData = Match_Metadata(MetadataList = self$Omics$Metadata)
      MetaDataColumn = as.vector(MetaData[MetadataFeature])[[1]]

      for(method in Methods){
        if(method == "SNF"){
          Distances = names(self$Multi_Omics$ClusterRes$SNF)
          for(distance in Distances){
            if(length(self$Multi_Omics$ClusterRes$SNF[[distance]]) < 2){
              next
            } else{
              ClusterRes = self$Multi_Omics$ClusterRes$SNF[[distance]]
              Title = paste0("Sankey plot SNF ", distance)
              Sankey_List_SNF = Data_Manipulation_Sankey_Clusters(ClusterResList = ClusterRes,
                                                                  MetadataColumn = MetaDataColumn)
              Sankey_SNF = Sankey_Plot(SankeyList = Sankey_List_SNF,
                                       Title = Title)
              self$Plots$Sankey_Clusters$SNF[[distance]] = Sankey_SNF
            }
          }
        } else if(method == "GAUDI"){
          GAUDIClustMethod = names(self$Multi_Omics$ClusterRes$GAUDI)
          for(clustmethod in GAUDIClustMethod){
            if(length(self$Multi_Omics$ClusterRes$GAUDI[[clustmethod]]) < 2){
              next
            } else{
              ClusterRes = self$Multi_Omics$ClusterRes$GAUDI[[clustmethod]]
              Title = paste0("Sankey plot GAUDI, clustered by ", clustmethod)
              Sankey_List_GAUDI = Data_Manipulation_Sankey_Clusters(ClusterResList = ClusterRes,
                                                                    MetadataColumn = MetaDataColumn)
              Sankey_GAUDI = Sankey_Plot(SankeyList = Sankey_List_GAUDI,
                                         Title = Title)
              self$Plots$Sankey_Clusters$GAUDI[[clustmethod]] = Sankey_GAUDI
            }
          }
        } else{
          if(length(self$Multi_Omics$ClusterRes[[method]]) < 2){
            next
          } else{
            ClusterRes = self$Multi_Omics$ClusterRes[[method]]
            Title = paste0("Sankey plot ", method)
            Sankey_List = Data_Manipulation_Sankey_Clusters(ClusterResList = ClusterRes,
                                                            MetadataColumn = MetaDataColumn)
            Sankey = Sankey_Plot(SankeyList = Sankey_List,
                                 Title = Title)
            self$Plots$Sankey_Clusters[[method]] = Sankey
          }
        }
      }
    },

    #Create Sankey plot comparing methods per cluster number split by SNF
    #' @description
        #' Interactive Sankey plot of sample cluster results between methods.
        #' Compares the cluster assignments of all multi-omics integration
        #' methods which have the same number of clusters. Can only be performed
        #' if multiple integration methods have calculated cluster assignments
        #' for the same cluster number.
    #' @param MetadataFeature String of the sample metadata feature for which the
    #' samples should be colored. The metadata feature need to be available
    #' in the metadata files for all omics types.
    #' @param Ensemble Whether to include the ensemble sample clustering results
    #' to the Sankey plot. Must run `run_Ensemble_Sample_MDS` or
    #' `run_Ensemble_Sample_CCA` first. Default is `TRUE`.
    #' @param SNFDistance Which SNF distance calculated cluster assignment
    #' should be included in the plot comparing different methods. To compare
    #' the cluster assignments between SNF distances, an SNF-only Sankey plot
    #' can be created, see parameter `AllPlots`. Must be one of:
    #' `"euclidean squared"` (default), `"euclidean"`, `"manhattan"`,
    #' `"minkowksi 0.25"`, `"minkowski 0.5"`, `"minkowski 3"` or `"minkowski 4"`.
    #' @param AllPlots Whether 3 plots should be created per cluster number:
    #'   * Selected: All methods + specified `SNFDistance` cluster assignments.
    #'   * SNF: All SNF distance cluster assignments.
    #'   * Full: All methods cluster assignments with all SNF distances.
    #' Default is `TRUE`. If `FALSE`, only the first option, with all methods
    #' and the specified `SNFDistance` will be created.
    #' @returns A list of plots. Sankey plot with the cluster assignments of
    #' all multi-omics integration methods and the specified `SNFDistance`,
    #' stored by cluster number in $Plots$Sankey_Methods$Selected. If `AllPlots`
    #' was `TRUE`, Sankey plot with the cluster assignments of all calculated
    #' SNF distances, stored by cluster number in $Plots$Sankey_Methods$SNF. If
    #' `AllPlots` was `TRUE`, Sankey plot with the cluster assignment of all
    #' multi-omics integration methods and all SNF distances, stored by cluster
    #' number in $Plots$Sankey_Methods$Full.
    plot_Sankey_Methods = function(MetadataFeature = NULL, #Column of metadata used for coloring samples
                                   Ensemble = TRUE, #Whether to add consensus clustering to Sankey
                                   SNFDistance = "euclidean squared", #SNF distance to add to rest plot: euclidean, euclidean squared, manhattan, minkowski 0.25, minkowski 0.5, minkowski 3, minkowski 4
                                   AllPlots = TRUE #whether to store all 3 plots (all methods + 1 SNF distance, all SNF distances, only non-SNF methods), or just all methods + 1 SNF
    ){

      #List of clustering methods used
      Methods = names(self$Multi_Omics$ClusterRes)

      Sankey_Split_SNF_Data = list()
      Sankey_Split_Rest_Data = list()
      Sankey_Split_Full_Data = list()

      #Create SNF and rest result lists
      for(method in Methods){
        if(method != "SNF" & method != "GAUDI"){
          ClustNum = names(self$Multi_Omics$ClusterRes[[method]])
          for(clust in ClustNum){
            Sankey_Split_Rest_Data[[clust]][[method]] = self$Multi_Omics$ClusterRes[[method]][[clust]]
            Sankey_Split_Full_Data[[clust]][[method]] = self$Multi_Omics$ClusterRes[[method]][[clust]]
          }
        }else if(method == "GAUDI"){
          ClustMethods = names(self$Multi_Omics$ClusterRes$GAUDI)
          for(clustmethod in ClustMethods){
            ClustNum = names(self$Multi_Omics$ClusterRes$GAUDI[[clustmethod]])
            for(clust in ClustNum){
              Sankey_Split_GAUDI_Name = paste0(method, "_", clustmethod)
              Sankey_Split_Rest_Data[[clust]][[Sankey_Split_GAUDI_Name]] = self$Multi_Omics$ClusterRes$GAUDI[[clustmethod]][[clust]]
              Sankey_Split_Full_Data[[clust]][[Sankey_Split_GAUDI_Name]] = self$Multi_Omics$ClusterRes$GAUDI[[clustmethod]][[clust]]
            }
          }
        }else if(method == "SNF"){
          Distances = names(self$Multi_Omics$ClusterRes$SNF)
          for(distance in Distances){
            ClustNum = names(self$Multi_Omics$ClusterRes$SNF[[distance]])
            for(clust in ClustNum){
              Sankey_Split_SNF_Name = paste0(method, "_", distance)
              Sankey_Split_SNF_Data[[clust]][[Sankey_Split_SNF_Name]] = self$Multi_Omics$ClusterRes$SNF[[distance]][[clust]]
              Sankey_Split_Full_Data[[clust]][[Sankey_Split_SNF_Name]] = self$Multi_Omics$ClusterRes$SNF[[distance]][[clust]]

            }
          }
        }
      }

      #Add user input SNF distance to Rest list
      Sankey_Split_Rest_ClustNum = names(Sankey_Split_Rest_Data)
      for(clustnum in Sankey_Split_Rest_ClustNum){
        if(SNFDistance == "euclidean"){
          SNF_Sankey_Name = "SNF_euclidean"
          Sankey_Split_Rest_Data[[clustnum]][[SNF_Sankey_Name]] = Sankey_Split_SNF_Data[[clustnum]][[SNF_Sankey_Name]]
        } else if(SNFDistance == "euclidean squared"){
          SNF_Sankey_Name = "SNF_euclidean squared"
          Sankey_Split_Rest_Data[[clustnum]][[SNF_Sankey_Name]] = Sankey_Split_SNF_Data[[clustnum]][[SNF_Sankey_Name]]
        } else if(SNFDistance == "manhattan"){
          SNF_Sankey_Name = "SNF_manhattan"
          Sankey_Split_Rest_Data[[clustnum]][[SNF_Sankey_Name]] = Sankey_Split_SNF_Data[[clustnum]][[SNF_Sankey_Name]]
        } else if(SNFDistance == "minkowski 0.25"){
          SNF_Sankey_Name = "SNF_minkowski 0.25"
          Sankey_Split_Rest_Data[[clustnum]][[SNF_Sankey_Name]] = Sankey_Split_SNF_Data[[clustnum]][[SNF_Sankey_Name]]
        } else if(SNFDistance == "minkowski 0.5"){
          SNF_Sankey_Name = "SNF_minkowksi 0.5"
          Sankey_Split_Rest_Data[[clustnum]][[SNF_Sankey_Name]] = Sankey_Split_SNF_Data[[clustnum]][[SNF_Sankey_Name]]
        } else if(SNFDistance == "minkowski 3"){
          SNF_Sankey_Name = "SNF_minkowski 3"
          Sankey_Split_Rest_Data[[clustnum]][[SNF_Sankey_Name]] = Sankey_Split_SNF_Data[[clustnum]][[SNF_Sankey_Name]]
        } else if(SNFDistance == "minkowski 4"){
          SNF_Sankey_Name = "SNF_minkowski 4"
          Sankey_Split_Rest_Data[[clustnum]][[SNF_Sankey_Name]] = Sankey_Split_SNF_Data[[clustnum]][[SNF_Sankey_Name]]
        }
      }

      #Add ensemble clustering to rest and full plots
      if(Ensemble == TRUE){
        #Rest ensemble
        MDSClustNum = names(self$Ensemble$Samples$MDS$ClusterRes)
        for(clustnum in MDSClustNum){
          Sankey_Split_Rest_Data[[clustnum]]$EnsembleMDS = self$Ensemble$Samples$MDS$ClusterRes[[clustnum]]
          Sankey_Split_Full_Data[[clustnum]]$EnsembleMDS = self$Ensemble$Samples$MDS$ClusterRes[[clustnum]]
        }
        CCAFactNum = names(self$Ensemble$Samples$CCA)
        for(factnum in CCAFactNum){
          CCAClustNum = names(self$Ensemble$Samples$CCA[[factnum]]$ClusterRes)
          for(clustnum in CCAClustNum){
            Sankey_Split_Rest_Data[[clustnum]]$EnsembleCCA = self$Ensemble$Samples$CCA[[factnum]]$ClusterRes[[clustnum]]
            Sankey_Split_Full_Data[[clustnum]]$EnsembleCCA = self$Ensemble$Samples$CCA[[factnum]]$ClusterRes[[clustnum]]
          }
        }
      }

      Sankey_Split_SNF_Data_List = list()
      Sankey_Split_Rest_Data_List = list()
      Sankey_Split_Full_Data_List = list()
      MetaData = Match_Metadata(MetadataList = self$Omics$Metadata)
      MetaDataColumn = as.vector(MetaData[MetadataFeature])[[1]]

      #Create Rest dataframes
      for(clust in Sankey_Split_Rest_ClustNum){
        Sankey_Split_Rest_Method = names(Sankey_Split_Rest_Data[[clust]])
        for(method in Sankey_Split_Rest_Method){
          Sankey_Split_Rest_Colname = paste0(method, "_", clust)
          Sankey_Split_Rest_Data_List[[clust]]["Rows"] = as.data.frame(row.names(Sankey_Split_Rest_Data[[clust]][[method]]))
          names(Sankey_Split_Rest_Data[[clust]][[method]])[names(Sankey_Split_Rest_Data[[clust]][[method]]) == "Cluster"] = Sankey_Split_Rest_Colname
          Sankey_Split_Rest_Data_List[[clust]] = cbind(Sankey_Split_Rest_Data_List[[clust]], Sankey_Split_Rest_Data[[clust]][[method]][Sankey_Split_Rest_Colname])
          row.names(Sankey_Split_Rest_Data_List[[clust]]) = NULL
          Sankey_Split_Rest_Data_List[[clust]] = tibble::column_to_rownames(Sankey_Split_Rest_Data_List[[clust]], "Rows")
        }
      }

      #Create Rest Sankey plot
      Sankey_Split_Rest_Data_List_red = purrr::discard(Sankey_Split_Rest_Data_List, ~any(ncol(.x) < 2))

      Sankey_Split_Rest_Data_List_clustnum = names(Sankey_Split_Rest_Data_List_red)
      for(clust in Sankey_Split_Rest_Data_List_clustnum){
        Sankey_Split_Rest_Single_Data = Sankey_Split_Rest_Data_List_red[[clust]]
        Sankey_Split_Rest_List = Data_Manipulation_Sankey_Methods(ClusterResData = Sankey_Split_Rest_Single_Data,
                                                                  MetadataColumn = MetaDataColumn)

        Sankey_Split_Rest_ncol = ncol(Sankey_Split_Rest_Single_Data)
        nclust = strsplit(clust, split = "_")[[1]][2]
        Title = paste0("Sankey plot comparing ", Sankey_Split_Rest_ncol, " methods with ", nclust, " clusters")
        Sankey_Split_Rest_Plot = Sankey_Plot(SankeyList = Sankey_Split_Rest_List,
                                             Title = Title)
        self$Plots$Sankey_Methods$Selected[[clust]] = Sankey_Split_Rest_Plot
      }

      #Create SNF and Full Sankey plots if AllPlots == TRUE
      if(AllPlots == TRUE){
        #Create SNF dataframes
        Sankey_Split_SNF_ClustNum = names(Sankey_Split_SNF_Data)
        for(clust in Sankey_Split_SNF_ClustNum){
          Sankey_Split_SNF_Method = names(Sankey_Split_SNF_Data[[clust]])
          for(method in Sankey_Split_SNF_Method){
            Sankey_Split_SNF_Colname = paste(method, "_", clust)
            Sankey_Split_SNF_Data_List[[clust]]["Rows"] = as.data.frame(row.names(Sankey_Split_SNF_Data[[clust]][[method]]))
            names(Sankey_Split_SNF_Data[[clust]][[method]])[names(Sankey_Split_SNF_Data[[clust]][[method]]) == "Cluster"] = Sankey_Split_SNF_Colname
            Sankey_Split_SNF_Data_List[[clust]] = cbind(Sankey_Split_SNF_Data_List[[clust]], Sankey_Split_SNF_Data[[clust]][[method]][Sankey_Split_SNF_Colname])
            row.names(Sankey_Split_SNF_Data_List[[clust]]) = NULL
            Sankey_Split_SNF_Data_List[[clust]] = tibble::column_to_rownames(Sankey_Split_SNF_Data_List[[clust]], "Rows")
          }
        }

        #Create SNF Sankey plot
        Sankey_Split_SNF_Data_List_red = purrr::discard(Sankey_Split_SNF_Data_List, ~any(ncol(.x) < 2))
        Sankey_Split_SNF_Data_List_clustnum = names(Sankey_Split_SNF_Data_List_red)
        for(clust in Sankey_Split_SNF_Data_List_clustnum){
          Sankey_Split_SNF_Single_Data = Sankey_Split_SNF_Data_List_red[[clust]]
          Sankey_Split_SNF_List = Data_Manipulation_Sankey_Methods(ClusterResData = Sankey_Split_SNF_Single_Data,
                                                                   MetadataColumn = MetaDataColumn)

          Sankey_Split_SNF_ncol = ncol(Sankey_Split_SNF_Single_Data)
          nclust = strsplit(clust, split = "_")[[1]][2]
          Title = paste0("Sankey plot comparing ", Sankey_Split_SNF_ncol, " SNF distances with ", nclust, " clusters")
          Sankey_Split_SNF_Plot = Sankey_Plot(SankeyList = Sankey_Split_SNF_List,
                                              Title = Title)
          self$Plots$Sankey_Methods$SNF[[clust]] = Sankey_Split_SNF_Plot
        }


        #Create Full dataframes
        Sankey_Split_Full_ClustNum = names(Sankey_Split_Full_Data)
        for(clust in Sankey_Split_Full_ClustNum){
          Sankey_Split_Full_Method = names(Sankey_Split_Full_Data[[clust]])
          for(method in Sankey_Split_Full_Method){
            Sankey_Split_Full_Colname = paste0(method, "_", clust)
            Sankey_Split_Full_Data_List[[clust]]["Rows"] = as.data.frame(row.names(Sankey_Split_Full_Data[[clust]][[method]]))
            names(Sankey_Split_Full_Data[[clust]][[method]])[names(Sankey_Split_Full_Data[[clust]][[method]]) == "Cluster"] = Sankey_Split_Full_Colname
            Sankey_Split_Full_Data_List[[clust]] = cbind(Sankey_Split_Full_Data_List[[clust]], Sankey_Split_Full_Data[[clust]][[method]][Sankey_Split_Full_Colname])
            row.names(Sankey_Split_Full_Data_List[[clust]]) = NULL
            Sankey_Split_Full_Data_List[[clust]] = tibble::column_to_rownames(Sankey_Split_Full_Data_List[[clust]], "Rows")
          }
        }

        #Create Full Sankey plot
        Sankey_Split_Full_Data_List_red = purrr::discard(Sankey_Split_Full_Data_List, ~any(ncol(.x) < 2))

        Sankey_Split_Full_Data_List_clustnum = names(Sankey_Split_Full_Data_List_red)
        for(clust in Sankey_Split_Full_Data_List_clustnum){
          Sankey_Split_Full_Single_Data = Sankey_Split_Full_Data_List_red[[clust]]
          Sankey_Split_Full_List = Data_Manipulation_Sankey_Methods(ClusterResData = Sankey_Split_Full_Single_Data,
                                                                    MetadataColumn = MetaDataColumn)

          Sankey_Split_Full_ncol = ncol(Sankey_Split_Full_Single_Data)
          nclust = strsplit(clust, split = "_")[[1]][2]
          Title = paste0("Sankey plot comparing ", Sankey_Split_Full_ncol, " methods with ", nclust, " clusters")
          Sankey_Split_Full_Plot = Sankey_Plot(SankeyList = Sankey_Split_Full_List,
                                               Title = Title)
          self$Plots$Sankey_Methods$Full[[clust]] = Sankey_Split_Full_Plot
        }
      }
    },

    #Top feature weights
    #' @description
        #' Plot the features with the highest absolute weights in each factor
        #' (component), per omics type and multi-omics integration method.
    #' @param nFeatures Number of features to plot for each factor. Default is
    #' 10.
    #' @param Scale Whether to scale the data from the different omics types, to
    #' make the weights more comparable between omics types. Default is `TRUE`.
    #' @returns A list of plots. Top features for each factor (component),
    #' stored by multi-omics integration method, number of factors and omics
    #' type in Plots$Top_Feature_Weights.
    plot_Top_Feature_Weights = function(nFeatures = 10,
                                        Scale = TRUE){
      Data = self$Multi_Omics$FeatureRes
      Feat_Weight_Methods = names(Data)

      for(method in Feat_Weight_Methods){
        if(method == "MoCluster" | method == "MCIA" | method == "jNMF" | method == "iNMF"){
          FactNum = names(Data[[method]])
          for(fact in FactNum){
            OmicsNames = names(Data[[method]][[fact]])
            for(omics in OmicsNames){
              names(Data[[method]][[fact]][[omics]]) = gsub(x = names(Data[[method]][[fact]][[omics]]), pattern = "V", replacement = "Factor")
              rownames(Data[[method]][[fact]][[omics]]) <- NULL
              Data[[method]][[fact]][[omics]] = tibble::column_to_rownames(Data[[method]][[fact]][[omics]], var = "Features")
              Data[[method]][[fact]][[omics]][] = sapply(Data[[method]][[fact]][[omics]], as.numeric)
              Plot = Top_Feature_Weights_All_Factors(FeatureRes = Data[[method]][[fact]][[omics]],
                                                     NumberFeatures = nFeatures,
                                                     Scale = Scale)
              self$Plots$Top_Feature_Weights[[method]][[fact]][[omics]] = Plot
            }
          }
        }else{
          FactNum = names(Data[[method]])
          for(fact in FactNum){
            OmicsNames = names(Data[[method]][[fact]])
            nfact = strsplit(fact, split = "_")[[1]][2]
            for(omics in OmicsNames){
              rownames(Data[[method]][[fact]][[omics]]) <- NULL
              Data[[method]][[fact]][[omics]] = tibble::column_to_rownames(Data[[method]][[fact]][[omics]], var = "Features")
              Data[[method]][[fact]][[omics]] = Data[[method]][[fact]][[omics]][1:nfact]
              Plot = Top_Feature_Weights_All_Factors(FeatureRes = Data[[method]][[fact]][[omics]],
                                                     NumberFeatures = nFeatures,
                                                     Scale = Scale)
              self$Plots$Top_Feature_Weights[[method]][[fact]][[omics]] = Plot
            }
          }
        }
      }
    },

    #Dendrogram
    #' @description
        #' Plot all previously created dendrograms with metadata features
        #' labeled in colored bars.
    #' @param Data Which previously created dendrogram to plot. Must be one of:
    #'   * `"Single"`: Single omics sample hierarchical clustering results. Will
    #'   plot the dendrograms for all omics types provided.
    #'   * `"Ensemble MDS"`: Ensemble MDS sample hierarchical clustering results.
    #'   * `"Ensemble Factor"`: Ensemble factor hierarchical clustering results,
    #'   as calculated with `run_Feature_Ensemble_HClust` with `Method`
    #'   `"Concatenation"`.
    #'   * `"Ensemble Feature Average"`: Ensemble feature hierarchical
    #'   clustering results, as calculated with `run_Feature_Ensemble_HClust`
    #'   with `Method` `"Average"`.
    #'   * `"Ensemble Feature Concatenation"`: Ensemble feature hierarchical
    #'   clustering results, as calculated with `run_Feature_Ensemble_HClust`
    #'   with `Method` `"Concatenation"`.
    #'   * `"Multi-Omics Feature"`: Feature hierarchical clustering results per
    #'   multi-omics integration method, as calculated with
    #'   `run_Multi_Omics_Feature_HClust`.
    #' @param Clusters Number of clusters to color the dendrogram branches for.
    #' If `NULL`, the branches will be black.
    #' @param MetadataFeatures Names of the sample metadata features to be
    #' included in colored bars below the dendrogram. A string vector with the
    #' names of the metadata features as in the metadata files should be
    #' provided. Can be used in `Data` = `"Single"` and `Data` =
    #' `"Ensemble MDS"`.
    #' @param LabelMetadata Name of the metadata feature for which the sample
    #' labels should be colored. Must be a sting with the name of the metadata
    #' feature as in the metadata files. Can be used in `Data` = `"Single"` and
    #' `Data` = `"Ensemble MDS"`.
    #' @param MetadataHeight Height of the colored bars with metadata
    #' information. Default is 0. If 0 puts the colored bars in the dendrogram,
    #' choose a negative integer. For `Data` = `"Ensemble MDS"`,
    #' `"Ensemble Feature Average"` or `"Ensemble Feature Concatenation"`, a
    #' single integer should be provided. For `Data` = `"Single"` or
    #' `"Multi-Omics Feature"`, a vector of integers with the same length as the
    #' number of omics datasets must be provided.
    #' @param LabelSize Size of the label beneath the dendrogram. Default is 1.
    #' @param BarLabelSize Size of the label of the colored bars. Default is 1.
    #' @returns A single or multiple dendrograms. If `Data` = `"Single"`, a
    #' sample dendrogram per omics type is returned, with colored bars with
    #' metadata features below each dendrogram. If `Data` = `"Ensemble MDS"`,
    #' a sample dendrogram is returned, with colored bars with metadata features
    #' below the dendrogram. If `Data` = `"Ensemble Factor"`, a factor
    #' dendrogram is returned. If `Data` =  `"Ensemble Feature Average"`, a
    #' feature dendrogram is returned, with a colored bar indicating the feature
    #' omics type below the dendrogram. If `Data` =
    #' `"Ensemble Feature Concatenation"`, a feature dendrogram is returned,
    #' with a colored bar indicating the feature omics type below the dendrogram.
    #' If `Data` = `"Multi-Omics Feature"`, a feature dendrogram per multi-omics
    #' integration method is returned, with colored bars indicating the feature
    #' omics type below each dendrogram.
    plot_Dendrogram = function(Data, #Single, Ensemble Sample, Ensemble Factor, Ensemble Feature Av, Ensemble Feature Concat
                               Clusters = NULL,
                               MetadataFeatures,
                               LabelMetadata,
                               MetadataHeight = 0,
                               LabelSize = 1,
                               BarLabelSize = 1
    ){
      if(Data == "Single"){
        if(is.null(self$Single_Omics$HClustTree)){
          stop("No single omics hierarchical clustering results, run run_Single_Omics_Hierarchical_Clustering first")
        } else{
          HCList = self$Single_Omics$HClustTree
          Omics = names(HCList)
          if(length(Omics) > length(MetadataHeight)){
            stop("Please provide MetadataHeight for all single omics types")
          }
          names(MetadataHeight) = Omics
          for(Omic in Omics){
            Dendro = stats::as.dendrogram(self$Single_Omics$HClustTree[[Omic]])
            Meta = self$Omics$Metadata[[Omic]]
            Title = paste0("Single omics hierarchical clustering of ", Omic)
            Colored_Meta = Dendrogram_Sample_Meta(MetadataFeatures = MetadataFeatures,
                                                  Metadata = Meta,
                                                  Dendrogram = Dendro,
                                                  LabelMetadata = LabelMetadata)
            if(is.null(Clusters)){
              Dendro |> dendextend::set("labels_colors", Colored_Meta$Label) |>
                dendextend::set("labels_cex", LabelSize) |>
                plot(main = Title)
              dendextend::colored_bars(colors = Colored_Meta$Vector,
                                       dend = Dendro,
                                       rowLabels = colnames(Colored_Meta$Vector),
                                       cex.rowLabels = BarLabelSize,
                                       y_shift = MetadataHeight[Omic])
            }else{
              Dendro |> dendextend::set("labels_colors", Colored_Meta$Label) |>
                dendextend::set("branches_k_color", k = Clusters) |>
                dendextend::set("labels_cex", LabelSize) |>
                plot(main = Title)
              dendextend::colored_bars(colors = Colored_Meta$Vector,
                                       dend = Dendro,
                                       rowLabels = colnames(Colored_Meta$Vector),
                                       cex.rowLabels = BarLabelSize,
                                       y_shift = MetadataHeight)
            }

          }
        }
      } else if(Data == "Ensemble MDS"){
        if(is.null(self$Plots$Dendrogram$Ensemble$MDS)){
          stop("No ensemble clustering results, run run_Ensemble_MDS first")
        } else{
          DendroList = self$Plots$Dendrogram$Ensemble$MDS
          Meta = Match_Metadata(self$Omics$Metadata)
          Clusterings = names(DendroList)
          for(Clust in Clusterings){
            Title = "Ensemble MDS sample clustering"
            Dendro = DendroList[[Clust]]
            Colored_Meta = Dendrogram_Sample_Meta(MetadataFeatures = MetadataFeatures,
                                                  Metadata = Meta,
                                                  Dendrogram = Dendro,
                                                  LabelMetadata = LabelMetadata)
            if(is.null(Clusters)){
              Dendro |> dendextend::set("labels_colors", Colored_Meta$Label) |>
                dendextend::set("labels_cex", LabelSize) |>
                plot(main = Title)
              dendextend::colored_bars(colors = Colored_Meta$Vector,
                                       dend = Dendro,
                                       rowLabels = colnames(Colored_Meta$Vector),
                                       cex.rowLabels = BarLabelSize,
                                       y_shift = MetadataHeight)
            }else{
              Dendro |> dendextend::set("labels_colors", Colored_Meta$Label) |>
                dendextend::set("branches_k_color", k = Clusters) |>
                dendextend::set("labels_cex", LabelSize) |>
                plot(main = Title)
              dendextend::colored_bars(colors = Colored_Meta$Vector,
                                       dend = Dendro,
                                       rowLabels = colnames(Colored_Meta$Vector),
                                       cex.rowLabels = BarLabelSize,
                                       y_shift = MetadataHeight)
            }
          }
        }
      } else if(Data == "Ensemble Factor"){
        if(is.null(self$Plots$Dendrogram$Factors$Concatenation)){
          stop("No ensemble factor results, run run_Feature_Ensemble_HClust with Method = Concatenation first")
        } else{
          Dendro = self$Plots$Dendrogram$Factors$Concatenation
          HClust = self$Ensemble$Factors$Concatenation$HClust
          Labels = Dendrogram_Label_Groups(Feature = F,
                                           HClust = HClust,
                                           Dendrogram = Dendro,
                                           LabelSize = F)
          Title = "Ensemble factor clustering"
          if(is.null(Clusters)){
            Dendro |>
              dendextend::set("labels_colors", Labels$Colors) |>
              dendextend::set("labels_cex", LabelSize) |>
              plot(main = Title)
            graphics::legend("topright", legend = levels(Labels$Legend_Levels), fill = Labels$Legend_Colors)
          }else{
            Dendro |>
              dendextend::set("labels_colors", Labels$Colors) |>
              dendextend::set("labels_cex", LabelSize) |>
              dendextend::set("branches_k_color", k = Clusters) |>
              plot(main = Title)
            graphics::legend("topright", legend = levels(Labels$Legend_Levels), fill = Labels$Legend_Colors)
          }
        }
      } else if(Data == "Ensemble Feature Average"){
        if(is.null(self$Plots$Dendrogram$Features$Average)){
          stop("No ensemble feature method average results, run run_Feature_Ensemble_HClust with Method = Average first")
        } else{
          Dendro = self$Plots$Dendrogram$Features$Average
          HClust = self$Ensemble$Features$Average$HClust
          Colored_Omics = Dendrogram_Omics_Bar(HClust = HClust,
                                               Dendrogram = Dendro)
          Title = "Ensemble feature clustering, method average"
          if(is.null(Clusters)){
            Dendro |>
              dendextend::set("labels", NA) |>
              plot(main = Title)
            dendextend::colored_bars(colors = Colored_Omics$Vector,
                                     dend = Dendro,
                                     rowLabels = "Omics",
                                     cex.rowLabels = BarLabelSize,
                                     y_shift = MetadataHeight)
            graphics::legend("topright", legend = levels(Colored_Omics$Legend_Levels), fill = Colored_Omics$Legend_Colors)
          }else{
            Dendro |>
              dendextend::set("labels", NA) |>
              dendextend::set("branches_k_color", k = Clusters) |>
              plot(main = Title)
            dendextend::colored_bars(colors = Colored_Omics$Vector,
                                     dend = Dendro,
                                     rowLabels = "Omics",
                                     cex.rowLabels = BarLabelSize,
                                     y_shift = MetadataHeight)
            graphics::legend("topright", legend = levels(Colored_Omics$Legend_Levels), fill = Colored_Omics$Legend_Colors)
          }
        }
      } else if(Data == "Ensemble Feature Concatenation"){
        if(is.null(self$Plots$Dendrogram$Features$Concatenation)){
          stop("No ensemble feature concatenation results, run run_Feature_Ensemble_HClust with Method = Concatenation first")
        } else{
          Dendro = self$Plots$Dendrogram$Features$Concatenation
          HClust = self$Ensemble$Features$Concatenation$HClust
          Colored_Omics = Dendrogram_Omics_Bar(HClust = HClust,
                                               Dendrogram = Dendro)
          Title = "Ensemble feature clustering, concatenation"
          if(is.null(Clusters)){
            Dendro |>
              dendextend::set("labels", NA) |>
              plot(main = Title)
            dendextend::colored_bars(colors = Colored_Omics$Vector,
                                     dend = Dendro,
                                     rowLabels = "Omics",
                                     cex.rowLabels = BarLabelSize,
                                     y_shift = MetadataHeight)
            graphics::legend("topright", legend = levels(Colored_Omics$Legend_Levels), fill = Colored_Omics$Legend_Colors)
          }else{
            Dendro |>
              dendextend::set("labels", NA) |>
              dendextend::set("branches_k_color", k = Clusters) |>
              plot(main = Title)
            dendextend::colored_bars(colors = Colored_Omics$Vector,
                                     dend = Dendro,
                                     rowLabels = "Omics",
                                     cex.rowLabels = BarLabelSize,
                                     y_shift = MetadataHeight)
            graphics::legend("topright", legend = levels(Colored_Omics$Legend_Levels), fill = Colored_Omics$Legend_Colors)
          }
        }
      } else if(Data == "Multi-Omics Feature"){
        if(is.null(self$Plots$Dendrogram$Features$Multi_Omics)){
          stop("No multi-omics feature hierarchical clustering results, run run_Multi_Omics_Feature_HClust first")
        } else{
          FeatMethods = names(self$Plots$Dendrogram$Features$Multi_Omics)
          for(method in FeatMethods){
            Dendro = self$Plots$Dendrogram$Features$Multi_Omics[[method]]
            HClust = self$Multi_Omics$Feature_HClust$HClust[[method]]
            Colored_Omics = Dendrogram_Omics_Bar(HClust = HClust,
                                                 Dendrogram = Dendro)
            Title = paste0("Multi-omics feature hierarchical clustering for ", method)

            if(is.null(Clusters)){
              Dendro |>
                dendextend::set("labels", NA) |>
                plot(main = Title)
              dendextend::colored_bars(colors = Colored_Omics$Vector,
                                       dend = Dendro,
                                       rowLabels = "Omics",
                                       cex.rowLabels = BarLabelSize,
                                       y_shift = MetadataHeight)
              graphics::legend("topright", legend = levels(Colored_Omics$Legend_Levels), fill = Colored_Omics$Legend_Colors)
            }else{
              Dendro |>
                dendextend::set("labels", NA) |>
                dendextend::set("branches_k_color", k = Clusters) |>
                plot(main = Title)
              dendextend::colored_bars(colors = Colored_Omics$Vector,
                                       dend = Dendro,
                                       rowLabels = "Omics",
                                       cex.rowLabels = BarLabelSize,
                                       y_shift = MetadataHeight)
              graphics::legend("topright", legend = levels(Colored_Omics$Legend_Levels), fill = Colored_Omics$Legend_Colors)
            }
          }
        }
      }else{
        stop("Unknown data type, please select one of the following: Single, Ensemble MDS, Ensemble Factor, Ensemble Feature Average, Ensemble Feature Concatenation, Multi-Omics Feature")
      }
    },

    #Overrepresentation analysis dot plot
    #' @description
        #' Plots a dotplot with the results from the overrepresentation analysis
        #' performed in `run_Over_Represenation`. Creates a dotplot with the
        #' overrepresented terms per omics type for each cluster of the ensemble
        #' feature clustering.
    #' @param Method Which ensemble feature clustering result is to be used.
    #' Must be one of" `"Average"` or `"Concatenation"`. Default is `"Average"`.
    #' @param Categories Maximum number of overrepresented GO terms to be
    #' displayed per cluster. Default is 10.
    #' @param FontSize Size of the font used in the plot. Default is 10.
    #' @returns A list of plots. A dotplot with the overrepresented GO terms per
    #' omics type for each cluster, stored by `Method`, omics type and cluster
    #' number in $Plots$ORA_Dot.
    plot_Enrich_Dot = function(Method = "Average", #Feature hierarchical clustering method, either Average or Concatenation
                               Categories = 10,
                               FontSize = 10){
      if(Method == "Average" | Method == "average" | Method == "Av" | Method == "av"){
        if(is.null(self$Ensemble$Features$Average$ORA)){
          stop("No feature clustering available for method average, please run run_Feature_Clustering with method Average")
        }else{
          EnrichData = self$Ensemble$Features$Average$ORA
        }
      }else if(Method == "Concatenation" | Method == "concatenation" | Method == "Concat" | Method == "concat" | Method == "Con" | Method == "con"){
        if(is.null(self$Ensemble$Features$Concatenation$ORA)){
          stop("No feature clustering available for concatenation, please run run_Feature_Clustering with method Concatenation")
        }else{
          EnrichData = self$Ensemble$Features$Concatenation$ORA
        }
      }else{
        stop("Unknown method of feature hierarchical clustering, please select either Average or Concatenation")
      }

      EnrichPlots = list()

      OmicsNames = names(EnrichData)
      for(omics in OmicsNames){
        EnrichClusters = names(EnrichData[[omics]])
        EnrichClusters = EnrichClusters[!EnrichClusters == "FeatureList"]
        for(cluster in EnrichClusters){
          ClusterEnrich = EnrichData[[omics]][[cluster]]$Enrichment
          if(is.null(ClusterEnrich)) next
          if(nrow(ClusterEnrich@result) == 0) next
          EnrichTitle = paste0("Over representation analysis of ", omics, " ", cluster)
          EnrichDotPlot = enrichplot::dotplot(object = ClusterEnrich,
                                              x = "GeneRatio",
                                              color = "p.adjust",
                                              showCategory = Categories,
                                              font.size = FontSize,
                                              title = EnrichTitle)
          EnrichPlots[[omics]][[cluster]] = EnrichDotPlot
        }
      }

      if(Method == "Average" | Method == "average" | Method == "Av" | Method == "av"){
        self$Plots$ORA_Dot$Average = EnrichPlots
      }else if(Method == "Concatenation" | Method == "concatenation" | Method == "Concat" | Method == "concat" | Method == "Con" | Method == "con"){
        self$Plots$ORA_Dot$Concatenation = EnrichPlots
      }
    }
  )
)
