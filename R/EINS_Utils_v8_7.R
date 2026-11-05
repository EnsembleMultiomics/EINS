#----------------Data processing----------------
#' Determine file delimiter
#' @description
#' Determine file delimiter of file in file path, used in EINS R6 class function
#' `add_Omics_file`.
#' @param path File path name for the data file.
#' @returns Vector with file delimiter.
#' @keywords internal
#' @references Olivier-Jimenez, D., Derks, R.J.E., Harari, O., Cruchaga, C.,
#' Ali, M., Ori, A., Di Fraia, D., Cabukusta, B., Henrie, A., Giera, M., and
#' Mohammed, Y. (2025). iSODA: A Comprehensive Tool for Integrative Omics Data
#' Analysis in Single- and Multi-Omics Experiments. Analytical Chemistry,
#' Vol. 97, Issue 5, Pages 2689-2697.
find_delim = function(path) {
  probe = paste(readLines(con = path, n = 10), collapse = "")
  sep = c("\t" = lengths(regmatches(probe, gregexpr("\t", probe))),
          "," = lengths(regmatches(probe, gregexpr(",", probe))),
          ";" = lengths(regmatches(probe, gregexpr(";", probe))))
  return(names(which.max(sep)))
}

#' Import data files from file path
#' @description
#' Import data files from file path, used in EINS R6 class function
#' `add_Omics_file`.
#' @param file_path File path name for the data file.
#' @param sep Separator used in the data file. Default is `NA`.
#' @param first_column_as_index Whether to use the first column of the data file
#' as the index. Default is `FALSE`.
#' @returns Data table from the file.
#' @keywords internal
#' @references Olivier-Jimenez, D., Derks, R.J.E., Harari, O., Cruchaga, C.,
#' Ali, M., Ori, A., Di Fraia, D., Cabukusta, B., Henrie, A., Giera, M., and
#' Mohammed, Y. (2025). iSODA: A Comprehensive Tool for Integrative Omics Data
#' Analysis in Single- and Multi-Omics Experiments. Analytical Chemistry,
#' Vol. 97, Issue 5, Pages 2689-2697.
soda_read_table = function(file_path,
                           sep = NA,
                           first_column_as_index = FALSE){

  if (is.na(sep)) {
    if (stringr::str_sub(file_path, -4, -1) == ".tsv") {
      sep = '\t'
    }
  }

  if (first_column_as_index) {
    index = 1
  } else {
    index = NULL
  }

  if (stringr::str_sub(file_path, -5, -1) == ".xlsx") {
    data_table = as.data.frame(readxl::read_xlsx(file_path))
  } else {
    if (is.na(sep)) {
      sep = find_delim(path = file_path)
      data_table = utils::read.csv(file_path,
                                   header = T,
                                   sep = sep,
                                   check.names = FALSE)
    } else {
      data_table = utils::read.csv(file_path,
                                   header = T,
                                   sep = sep,
                                   check.names = FALSE)
    }

  }

  if (!is.null(index)) {

    duplicates = duplicated(data_table[,index])
    if (sum(duplicates) > 0) {
      warning(paste0('Removed ', sum(duplicates), ' duplicated samples'))
      data_table = data_table[!duplicates,]
    }
    rownames(data_table) = data_table[,1]
    data_table[,1] = NULL
  }



  original_count = ncol(data_table)
  if (original_count > 1) {
    data_table = data_table[,!duplicated(colnames(data_table))]
    final_count = ncol(data_table)
    if(original_count != final_count) {
      warning(paste0('Removed ', original_count - final_count, ' duplicated columns'))
    }
  }
  return(data_table)
}

#Calculate feature proportion of feature total in sample
#' Calculate the feature proportion of the total feature amount per sample
#' @description
#' Calculate the feature proportion of the total feature amount per sample,
#' used in EINS R6 class function `run_Preprocessing`, if `FunctionOrder`
#' includes `"FeatureProp"`.
#' @param OmicSampleData Single omics dataset, with features in rows.
#' @returns Single omics matrix with proportion of each feature of feature total
#' in the sample calculated.
#' @keywords internal
Feature_Proportion = function(OmicSampleData #single omics dataset (samples in columns)
){
  Prop_Data = OmicSampleData/colSums(OmicSampleData, na.rm = TRUE)[col(OmicSampleData)]
  return(Prop_Data)
}

#filter out features with higher NA proportion than accepted
#' Filter features with NA proportions above user-defined threshold
#' @description
#' Filter features based on their proportion of NAs per feature using a
#' user-defined threshold. Used in EINS R6 class function `run_Preprocessing`,
#' if `FunctionOrder` includes `"FilterCoverage"`.
#' @param OmicSampleData Single omics dataset, with features in rows.
#' @param Coverage Numerical threshold for feature coverage. Default is 0.8
#' (at least 0.8 of the samples are not NA per feature).
#' @returns Single omics matrix with features with more NAs than the
#' threshold allows for filtered out.
#' @keywords internal
filter_Coverage = function(OmicSampleData, #single omics dataset samples only (samples in columns)
                           Coverage = 0.8 #proportion of samples with measurement to retain feature
){
  Data = OmicSampleData
  RemovedRows = c()
  #Calculate NA percentage per feature (in rows)
  NAPercent = rowMeans(is.na(Data))
  RemovedRows = names(which(NAPercent > (1 - Coverage)))

  SampleData = as.data.frame(Data)
  SampleDataRows = tibble::rownames_to_column(SampleData, var = "rowname")
  SampleDataRowsRem = SampleDataRows[!(SampleDataRows$rowname %in% RemovedRows), ]
  rownames(SampleDataRowsRem) = NULL
  SampleDataRowsRem = tibble::column_to_rownames(SampleDataRowsRem, var = "rowname")
  SampleDataRowsRem = as.matrix(SampleDataRowsRem)

  return(SampleDataRowsRem)
}

#filter out features where more than accepted threshold are below the blank mean * multiplier
#' Filter the features based on the measurement in blank samples
#' @description
#' Filter the features based on the measurement in blank samples. Multiply the
#' mean feature measurement from blank samples by a user-defined multiplier, and
#' filter out the features in which fewer features than the user-defined
#' threshold exceed the multiplied blank mean. Only possible if blank samples
#' are available. Used in EINS R6 class function `run_Preprocessing`, if
#' `FunctionOrder` includes `"FilterBlankMean"`.
#' @param OmicSampleData Single omics dataset, with features in rows.
#' @param OmicBlankData Single omics blank samples dataset, with features in
#' rows.
#' @param BlankMeanMultiplier Numerical by which the blank feature mean is
#' multiplied. Default is 2.
#' @param SampleThreshold Numerical threshold for blank mean filtering. Default
#' is 0.8 (at least 0.8 of the samples are larger than the blank feature mean *
#' `BlankMeanMultiplier`).
#' @returns Single omics matrix with features with more samples below the
#' multiplied blank mean than the threshold allows for filtered out.
#' @keywords internal
filter_Blank_Mean = function(OmicSampleData, #single omics dataset samples only (samples in columns)
                             OmicBlankData, #single omics dataset blanks only (samples in columns)
                             BlankMeanMultiplier = 2, #Multiplier of the blank mean to use as threshold
                             SampleThreshold = 0.8 #proportion of samples higher than Multiplier*blank mean to retain feature
){
  BlankMeans = rowMeans(OmicBlankData, na.rm = TRUE)
  RemovedRows = c()
  n_Samples = ncol(OmicSampleData)
  for(row in rownames(OmicSampleData)){
    BlankThreshold = BlankMeanMultiplier * BlankMeans[row]
    AboveThresholdCount = sum(OmicSampleData[row, ] > BlankThreshold, na.rm = TRUE)
    if((AboveThresholdCount/n_Samples) < SampleThreshold){
      RemovedRows = c(RemovedRows, row)
    }
  }

  SampleData = as.data.frame(OmicSampleData)
  SampleDataRows = tibble::rownames_to_column(SampleData, var = "rowname")
  SampleDataRowsRem = SampleDataRows[!(SampleDataRows$rowname %in% RemovedRows), ]
  rownames(SampleDataRowsRem) = NULL
  SampleDataRowsRem = tibble::column_to_rownames(SampleDataRowsRem, var = "rowname")
  SampleDataRowsRem = as.matrix(SampleDataRowsRem)

  return(SampleDataRowsRem)
}

#Z-score normalization
#' Feature Z-score normalization
#' @description
#' Perform Z-score normalization on the features. Used in EINS R6 class function
#' `run-Preprocessing`, if `FunctionOrder` includes `"Normalization"`.
#' @param OmicSampleData Single omics dataset, with features in rows.
#' @returns Single omics matrix with features normalized with Z-score
#' normalization.
#' @keywords internal
normalize_ZScore = function(OmicSampleData #single omics dataset samples only (samples in columns)
){
  RowMeans = rowMeans(OmicSampleData, na.rm = TRUE)
  RowSD = apply(OmicSampleData, 1, sd, na.rm = TRUE)
  MeanSubData = sweep(OmicSampleData, 1, RowMeans, "-")
  NormData = sweep(MeanSubData, 1, RowSD, "/")

  return(NormData)
}

#NA imputation
#' Missing value (NA) imputation
#' @description
#' NA imputation or replacement with one of the following methods: k-nearest
#' neighbor imputation, missForest imputation, zero replacement or NA omittance.
#' Used in EINS R6 class function `run_Preprocessing`, if `FunctionOrder`
#' includes `"NAImpute"`.
#' @param OmicSampleData Single omics dataset, with features in rows.
#' @param Method String with method of NA imputation. Must be one of: `"KNN"`,
#' `"missForest"`, `"zero replacement"`, or `"omit NA"`.
#' @param K_Neighbors Numerical, number of k nearest neighbors. Only required
#' when `Method` = `"KNN"`.
#' @returns Single omics matrix with NA values imputed, replaced, or features
#' omitted.
#' @keywords internal
#' @references   * Hastie T, Tibshirani R, Narasimhan B, Chu G (2025). impute:
#'   impute: Imputation for microarray data. doi:10.18129/B9.bioc.impute,
#'   R package version 1.82.0
#'   * Stekhoven DJ (2022). missForest: Nonparametric Missing Value Imputation
#'   using Random Forest. R package version 1.5.
impute_NA = function(OmicSampleData,
                     Method = NULL,
                     K_Neighbors = 10){
  #KNN
  if(Method == "KNN" | Method == "knn"){
    if(anyNA(OmicSampleData, recursive = TRUE) == TRUE){
      NAImputedData = impute::impute.knn(OmicSampleData, k = K_Neighbors)$data
    }else{
      NAImputedData = OmicSampleData
    }
    #missForest
  } else if(Method == "missForest" | Method == "missforest"){
    if(anyNA(OmicSampleData, recursive = TRUE) == TRUE){
      NAImputedData = missForest::missForest(OmicSampleData)$ximp
    }else{
      NAImputedData = OmicSampleData
    }
    #Zero replacement
  } else if(Method == "zero" | Method == "Zero" | Method == "Zero Replacement" | Method == "zero replacement"){
    if(anyNA(OmicSampleData, recursive = TRUE) == TRUE){
      NAImputedData = OmicSampleData |> replace(is.na(.), 0)
    }else{
      NAImputedData = OmicSampleData
    }
    #Omit NA
  } else if(Method == "omit" | Method == "Omit" | Method == "omit na" | Method == "Omit NA"){
    if(anyNA(OmicSampleData, recursive = TRUE) == TRUE){
      NAImputedData = as.data.frame(na.omit(OmicSampleData))
      NAImputedData = as.matrix(NAImputedData)
      if(nrow(OmicSampleData) != nrow(NAImputedData)){
        message(paste0((nrow(OmicSampleData) - nrow(NAImputedData)), " features with NA values are removed"))
      }
    }else{
      NAImputedData = OmicSampleData
    }
  } else{
    if(anyNA(OmicSampleData, recursive = TRUE) == TRUE){
      base::stop("Missing values in Omics data, select one of following methods to deal with NA: KNN, missForest, Zero Replacement, Omit NA")
    }else{
      NAImputedData = OmicSampleData
    }
  }
  return(NAImputedData)
}

#Batch effect PCA plot
#' PCA plot of sample batch effect
#' @description
#' Creates a PCA plot of the samples, colored by batch to examine the batch
#' effect. Only possible if batch information is available for the samples. Used
#' in EINS R6 class function `run_Preprocessing`, if `FunctionOrder` includes
#' `"BatchCorrection"`.
#' @param Data Single omics dataset, with features in rows.
#' @param Metadata Single omics metadata file, with batch information.
#' @param BatchColumnName String of the column name which contains batch
#' information in the metadata file.
#' @param XAxis Numerical to indicate which principle component will be used for
#' the X-axis in the PCA plot. Default is 1.
#' @param YAxis Numerical to indicate which principle component will be used for
#' the Y-acis in the PCA plot. Default is 2.
#' @param Title String, title for the batch effect PCA plot.
#' @returns Interactive scatterplot with samples colored by batch.
#' @keywords internal
PCA_Plot_Batch = function(Data,
                          Metadata,
                          BatchColumnName,
                          XAxis = 1,
                          YAxis = 2,
                          Title){
  #Extract PCA results
  PCA_rot = as.data.frame(Data$x)

  #Extract PC standard deviations
  PCA_sdev = Data$sdev
  #Calculate PC variance explained
  PCA_var_ex = Data$sdev^2/sum(Data$sdev^2)
  #Save PC variance explained for XAxis
  PCA_var_ex_per_X = round(PCA_var_ex[XAxis] * 100, 2)
  #Save PC variance explained for YAxis
  PCA_var_ex_per_Y = round(PCA_var_ex[YAxis] * 100, 2)
  #Extract the Batch information from the metadata file
  Metadata[BatchColumnName] = as.factor(Metadata[,BatchColumnName])

  Plot = plotly::plot_ly(data = PCA_rot,
                         x = PCA_rot[,XAxis],
                         y = PCA_rot[,YAxis],
                         color = Metadata[,BatchColumnName],
                         colors = "Set1")

  Plot = plotly::add_trace(p = Plot,
                           type = "scatter",
                           mode = "markers",
                           text = ~paste("Sample: ", rownames(PCA_rot))
  )
  Plot = plotly::layout(p = Plot,
                        title = list(text = Title),
                        xaxis = list(title = list(text = paste0(names(PCA_rot[XAxis]), " (", PCA_var_ex_per_X, "%)"),
                                                  font = list(size = 15))),
                        yaxis = list(title = list(text = paste0(names(PCA_rot[YAxis]), " (", PCA_var_ex_per_Y, "%)"),
                                                  font = list(size = 15))),
                        legend = list(title = list(text = "Batch"))
  )
  return(Plot)
}

#Batch effect correction function
#' Batch effect correction using ComBat from the sva package.
#' @description
#' Batch effect correction using ComBat from the sva package. Only possible if
#' batch information is available for the samples. Used in EINS R6 class
#' function `run_Preprocessing`, if `FunctionOrder` includes `"BatchCorrection"`.
#' @param OmicSampleData Single omics dataset, with features in rows.
#' @param Metadata single omics metadata file, with batch information.
#' @param BatchColumn String of the column name which contains batch information
#' in the metadata file.
#' @returns Single omics matrix with sample batch effect corrected.
#' @keywords internal
#' @references Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Zhang Y,
#' Storey JD, Torres LC (2025). sva: Surrogate Variable Analysis.
#' doi:10.18129/B9.bioc.sva, R package version 3.56.0
Batch_Effect_Correction = function(OmicSampleData,
                                   Metadata,
                                   BatchColumn
){
  Metadata = as.data.frame(Metadata)

  #create the model matrix
  modmat = stats::model.matrix(~1, data = Metadata)

  #create batch vector
  batches = Metadata[,BatchColumn]

  #correct for batch effect
  batchcor = sva::ComBat(dat = OmicSampleData,
                         batch = batches,
                         mod = modmat,
                         par.prior = TRUE,
                         prior.plots = FALSE,
                         mean.only = FALSE,
                         ref.batch = NULL)
  return(batchcor)
}

#Filter features based on pairwise correlation
#' Filter out correlated features
#' @description
#' Filter out correlated features from the omics dataset. Calculates pairwise
#' correlations between features, and filter out those features with correlation
#' above a user-defined threshold. Either one feature is kept as the
#' representative, with the filtered out features in a list by representative,
#' or the mean measurements from the correlated features is calculated and kept.
#' Used in EINS R6 class function `run_Preprocessing`, if `FunctionOrder`
#' includes `"FilterCorrelation"`.
#' @param OmicSampleData Single omics dataset, with features in rows.
#' @param CutOff Numerical threshold for correlation. Default is 0.9 (if
#' pairwise correlation between 2 features is 0.9 or larger, they are filtered).
#' @param CorMean Whether to calculate the mean measurements for the correlated
#' features. Default is `FALSE`. If `TRUE`, the mean of all feature measurements
#' for a group of correlated features is included in the dataset.
#' @returns Single omics matrix, with correlated features filtered out
#' @keywords internal
CorFilter = function(OmicSampleData , #omics matrix, features in rows
                     CutOff = 0.9, #correlation threshold
                     CorMean = FALSE #calculate sample mean of correlated group or take representative of group
){
  tData = t(OmicSampleData ) #features in columns
  if(anyNA(tData, recursive = TRUE) == TRUE){
    stop("The data contains missing values, put NAImpute before FilterCor in FunctionOrder")
  }

  #Feature correlation
  CorData = stats::cor(tData)
  CorData[lower.tri(CorData, diag = TRUE)] = NA
  CorData = replace(CorData, CorData >= abs(CutOff), 1)
  CorData = replace(CorData, CorData < abs(CutOff), 0)

  FeatList = list()
  for(i in 1:nrow(CorData)){
    if(i > nrow(CorData)){
      next
    }
    CorFeat = which(CorData[i,] == 1)
    if(length(CorFeat) != 0){
      FeatName = names(CorFeat)
      RowName = rownames(CorData)[i]
      FeatList[[RowName]] = c(RowName, FeatName)
      CorData = CorData[!rownames(CorData) %in% FeatName, !colnames(CorData) %in% FeatName]
    }
  }

  if(CorMean == TRUE){
    for(i in 1:length(FeatList)){
      Feats <- FeatList[[i]]
      MeanName <- paste0("MeanCor", i)
      GroupMean <- rowMeans(tData[, Feats])
      tData = tData[, !colnames(tData) %in% Feats]
      tData = cbind(tData, GroupMean)
      colnames(tData)[colnames(tData) == "GroupMean"] = MeanName
      names(FeatList)[i] <- MeanName
    }
  }else if(CorMean == FALSE){
    NoCorRow <- rownames(CorData)
    tData <- tData[, colnames(tData) %in% NoCorRow]
  }

  CorMeanData = t(tData) #samples back in column
  CorRes = list(RedMat = CorMeanData,
                CorGroups = FeatList)
  return(CorRes)
}

#Optimal cluster number with NbClust
#' Calculate optimal cluster number of integrated multi-omics data with NbClust.
#' @description
#' Calculate the optimal cluster number for integrated data using NbClust.
#' Performs data integration with seven methods that provide sample by
#' factor/component matrices: MoCluster, MCIA, jNMF, iNMF, LRAcluster, MOFA and
#' GAUDI. Integration methods are performed with default parameters. A range of
#' factors/components can be provided, for multiple integration results to be
#' used. NbClust is then performed on these matrices. Multiple NbClust indexes
#' are calculated per matrix. All results from NbClust indexes for all matrices
#' are combined and the mean, mode and median optimal cluster number are
#' calculated. The user can then determine the preferred cluster number. Used in
#' EINS R6 class function `run_Optimal_Cluster_Number`.
#' @param Data List of omics matrices, with samples in columns.
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
#' @keywords internal
#' @references   * Charrad M, Ghazzali N, Boiteau V, Niknafs A (2014). “NbClust:
#'   An R Package for Determining the Relevant Number of Clusters in a Data Set.”
#'   Journal of Statistical Software, 61(6), 1–36.
#'   * Meng C (2025). mogsa: Multiple omics data integrative clustering and gene
#'   set analysis. doi:10.18129/B9.bioc.mogsa, R package version 1.42.0
#'   * Meng C, Kuster B, Culhane A, Gholami AM (2013). “A multivariate approach
#'   to the integration of multi-omics datasets.” BMC Bioinformatics.
#'   * Tsuyuzaki, K., and Nikaido, I. (2024). nnTensor: Non-Negative Tensor
#'   Decomposition. R package version 1.3.0
#'   * Chalise, P., Raghavan, R., and Fridley, B. (2025). IntNMF: Integrative
#'   Clustering of Multiple Genomic Dataset. R package version 1.3.0
#'   * Lu, X., Meng, J., Zhou, Y., Jiang, L., and Yan, F. (2020). MOVICS: an R
#'   package for multi-omics integration and visualization in cancer subtyping.
#'   Bioinformatics, btaa1018.
#'   * Wu D, Wang D, Zhang MQ, Gu J (2015). Fast dimension reduction and
#'   integrative clustering of multi-omics data using low-rank approximation:
#'   application to cancer molecular classification. BMC Genomics, 16(1):1022.
#'   * Argelaguet R, Velten B, Arnol D, Dietrich S, Zenz T, Marioni JC, Buettner
#'   F, Huber W, Stegle O (2018). “Multi‐Omics Factor Analysis—a framework for
#'   unsupervised integration of multi‐omics data sets.” Molecular Systems
#'   Biology, 14.
#'   * Castellano-Escuder P, Zachman DK, Han K, Hirschey MD. GAUDI:
#'   interpretable multi-omics integration with UMAP embeddings and
#'   density-based clustering. Nat Commun. 2025 Jul 1;16(1):5771.
Optimal_Cluster_Number = function(Data = NULL, #list of omics data
                                  Components = NULL, #single value or range of possible components
                                  NbClustDistance = "euclidean",
                                  NbClustMinClust = 2,
                                  NbClustMaxClust = 12,
                                  NbClustMethod = "kmeans",
                                  NbClustIndex = "all",
                                  NbClustAlphaBeale = 0.1
){
  OptimalClusterRes = list()
  NbClustRes = list()
  NMFData = list()
  TData = list()
  #create NMF data
  NMFData = lapply(Data, function(dat){
    if(!all(dat >= 0)){
      dat = pmax(dat + abs(min(dat)), 0)
    }
    dat = dat/max(dat)
    dat = base::t(dat)
  })

  #create transposed data
  TData = lapply(Data, t)

  for(i in Components){
    ListName = paste0("Components_", i)
    #run MoCluster
    moaClust = mogsa::mbpca(x = Data,
                            ncomp = i,
                            k = "all",
                            method = "globalScore",
                            option = "lambda1",
                            center = TRUE,
                            scale = FALSE,
                            moa = TRUE,
                            svd.solver = "fast",
                            maxiter = 1000,
                            verbose = FALSE)
    MoClustFactScore = moaClust@fac.scr
    #NbClust on MoCluster
    NbClustMoCluster = NbClust::NbClust(data = MoClustFactScore,
                                        distance = NbClustDistance,
                                        min.nc = NbClustMinClust,
                                        max.nc = NbClustMaxClust,
                                        method = NbClustMethod,
                                        index = NbClustIndex,
                                        alphaBeale = NbClustAlphaBeale)
    NbClustRes$MoCluster[[ListName]] = NbClustMoCluster$Best.nc[1,]

    #run MCIA
    MCIAClust = omicade4::mcia(df.list = Data,
                               cia.nf = i,
                               cia.scan = FALSE,
                               nsc = TRUE,
                               svd = TRUE)
    MCIAClustFactScore = MCIAClust$mcoa$SynVar
    #NbClust on MCIA
    NbClustMCIA = NbClust::NbClust(data = MCIAClustFactScore,
                                   distance = NbClustDistance,
                                   min.nc = NbClustMinClust,
                                   max.nc = NbClustMaxClust,
                                   method = NbClustMethod,
                                   index = NbClustIndex,
                                   alphaBeale = NbClustAlphaBeale)
    NbClustRes$MCIA[[ListName]] = NbClustMCIA$Best.nc[1,]

    #run jNMF
    jNMFClust = nnTensor::jNMF(X = NMFData,
                               pseudocount = .Machine$double.eps,
                               J = i,
                               algorithm = "KL",
                               num.iter = 100)
    jNMFClustFactScore = jNMFClust$W
    #NbClust on jNMF
    NbClustjNMF = NbClust::NbClust(data = jNMFClustFactScore,
                                   distance = NbClustDistance,
                                   min.nc = NbClustMinClust,
                                   max.nc = NbClustMaxClust,
                                   method = NbClustMethod,
                                   index = NbClustIndex,
                                   alphaBeale = NbClustAlphaBeale)
    NbClustRes$jNMF[[ListName]] = NbClustjNMF$Best.nc[1,]

    #run iNMF
    iNMFClust = IntNMF::nmf.mnnals(dat = NMFData,
                                   k = i,
                                   maxiter = 200,
                                   st.count = 20,
                                   n.ini = 30)
    iNMFClustFactScore = iNMFClust$W
    #NbClust on iNMF
    NbClustiNMF = NbClust::NbClust(data = iNMFClustFactScore,
                                   distance = NbClustDistance,
                                   min.nc = NbClustMinClust,
                                   max.nc = NbClustMaxClust,
                                   method = NbClustMethod,
                                   index = NbClustIndex,
                                   alphaBeale = NbClustAlphaBeale)
    NbClustRes$iNMF[[ListName]] = NbClustiNMF$Best.nc[1,]

    #run LRAcluster
    LRAClust = getLRAcluster_MOVICS(data = Data,
                                    N.clust = i,
                                    type = rep("gaussian", length(Data)),
                                    clusterAlg = "complete")
    LRAClustFactScore = t(LRAClust$fit$coordinate)
    #NbClust on LRAcluster
    NbClustLRA = NbClust::NbClust(data = LRAClustFactScore,
                                  distance = NbClustDistance,
                                  min.nc = NbClustMinClust,
                                  max.nc = NbClustMaxClust,
                                  method = NbClustMethod,
                                  index = NbClustIndex,
                                  alphaBeale = NbClustAlphaBeale)
    NbClustRes$LRA[[ListName]] = NbClustLRA$Best.nc[1,]

    #run GAUDI
    Single_Omics_UMAP_Clust = list()
    for(j in 1:length(TData)){
      Single_Omics_UMAP_Clust[[i]] = uwot::umap(TData[[j]],
                                                n_neighbors = 10,
                                                n_components = i,
                                                metric = "euclidean")
    }
    SO_UMAP_Concat_Clust = dplyr::bind_cols(Single_Omics_UMAP_Clust, .name_repair = "unique_quiet")
    UMAPClustFactScore = uwot::umap(SO_UMAP_Concat_Clust,
                                    n_neighbors = 10,
                                    n_components = i,
                                    metric = "euclidean")
    #NbClust on GAUDI
    NbClustUMAP = NbClust::NbClust(data = UMAPClustFactScore,
                                   distance = NbClustDistance,
                                   min.nc = NbClustMinClust,
                                   max.nc = NbClustMaxClust,
                                   method = NbClustMethod,
                                   index = NbClustIndex,
                                   alphaBeale = NbClustAlphaBeale)
    NbClustRes$UMAP[[ListName]] = NbClustUMAP$Best.nc[1,]

    #run MOFA
    MOFAClustObj = MOFA2::create_mofa(data = Data,
                                      groups = NULL)

    MOFADataOpts = MOFA2::get_default_data_options(MOFAClustObj)
    MOFADataOpts$scale_views = FALSE

    MOFAModelOpts = MOFA2::get_default_model_options(MOFAClustObj)
    MOFAModelOpts$likelihoods = rep("gaussian", length(Data))
    MOFAModelOpts$num_factors = i
    MOFAModelOpts$spikeslab_factors = FALSE
    MOFAModelOpts$spikeslab_weights = FALSE
    MOFAModelOpts$ard_factors = FALSE
    MOFAModelOpts$ard_weights = TRUE

    MOFATrainOpts = MOFA2::get_default_training_options(MOFAClustObj)
    MOFATrainOpts$maxiter = 1000
    MOFATrainOpts$convergence_mode = "fast"
    MOFATrainOpts$startELBO = 1
    MOFATrainOpts$freqELBO = 1
    MOFATrainOpts$stochastic = FALSE

    MOFAClustObj = MOFA2::prepare_mofa(object = MOFAClustObj,
                                       data_options = MOFADataOpts,
                                       model_options = MOFAModelOpts,
                                       training_options = MOFATrainOpts)

    MOFAClustModel = MOFA2::run_mofa(object = MOFAClustObj,
                                     outfile = NULL,
                                     use_basilisk = TRUE,
                                     save_data = FALSE)

    MOFAClustFactScore = MOFAClustModel@expectations$Z$group1
    #NbClust on MOFA

    NbClustMOFA = NbClust::NbClust(data = MOFAClustFactScore,
                                   distance = NbClustDistance,
                                   min.nc = NbClustMinClust,
                                   max.nc = NbClustMaxClust,
                                   method = NbClustMethod,
                                   index = NbClustIndex,
                                   alphaBeale = NbClustAlphaBeale)
    NbClustRes$MOFA[[ListName]] = NbClustMOFA$Best.nc[1,]
  }

  #per method optimal cluster numbers by mean, median, most common

  OptimalClusterRes$NbIndexResults = NbClustRes
  NbClustResComb = list()
  NbMethods = names(NbClustRes)
  for(method in NbMethods){
    NbClustResComb[[method]] = unlist(NbClustRes[[method]])
    OptimalClusterRes$NbIndexMostCommon[[method]] = DescTools::Mode(unlist(NbClustRes[[method]]))[1]
    OptimalClusterRes$NbIndexMean[[method]] = round(mean(unlist(NbClustRes[[method]])), digits = 0)
    OptimalClusterRes$NbIndexMedian[[method]] = stats::median(unlist(NbClustRes[[method]]))
  }

  #select optimal cluster number from different methods by mean, median, most common
  OptimalClusterRes$Mean = round(mean(unlist(NbClustResComb)), digits = 0)
  OptimalClusterRes$Median = stats::median(unlist(NbClustResComb))
  OptimalClusterRes$MostCommon = DescTools::Mode(unlist(NbClustResComb))[1]
  return(OptimalClusterRes)
}

#----------------Multi-omics integration----------------

#MOVICS LRAcluster data management function
#' Wrapper function for LRAcluster, as used in getLRAcluster from the MOVICS
#' package.
#' @description
#' Wrapper function for LRAcluster, as used in getLRAcluster from the MOVICS
#' package. Used in the EINS R6 class function `run_LRAcluster`.
#' @param data List of omics matrices, with samples in columns.
#' @param N.clust Number of sample clusters to be created.
#' @param type String indicating the data type in the list of matrices. Must be
#' one of: `"gaussian"` (default), `"binomial"`, or `"possion"`.
#' @param clusterAlg Agglomeration method to be used for the hierarchical
#' clustering. Must be one of: `"ward.D"` (default), `"ward.D2"`, `"single"`,
#' `"complete"`, `"average"`, `"mcquitty"`, `"median"` or `"centroid"`.
#' @returns A list of results. $fit, the object return by `LRAcluster_MOVICS`.
#' $clust.res, a dataframe with the sample cluster results. $clust.dend, a
#' dendrogram of the sample clustering. $mo.method, a string value indicating
#' the method used for multi-omics integrative clustering
#' @keywords internal
#' @references Lu, X., Meng, J., Zhou, Y., Jiang, L., and Yan, F. (2020).
#' MOVICS: an R package for multi-omics integration and visualization in cancer
#' subtyping. Bioinformatics, btaa1018.
getLRAcluster_MOVICS <- function(data       = NULL,
                                 N.clust    = NULL,
                                 type       = rep("gaussian", length(data)),
                                 clusterAlg = "ward.D"){

  # check data
  n_dat <- length(data)
  if(n_dat > 6){
    stop('LRAcluster can support up to 6 datasets.')
  }
  if(n_dat < 2){
    stop('LRAcluster needs at least 2 omics data.')
  }

  data <- lapply(data, as.matrix)

  if(is.element("binomial",type)) {
    bindex <- which(type == "binomial")
    for (i in bindex) {
      a <- which(rowSums(data[[i]]) == 0)
      b <- which(rowSums(data[[i]]) == ncol(data[[i]]))
      if(length(a) > 0) {
        data[[i]] <- data[[i]][which(rowSums(data[[i]]) != 0),] # remove all zero
      }

      if(length(b) > 0) {
        data[[i]] <- data[[i]][which(rowSums(data[[i]]) != ncol(data[[i]])),] # remove all one
      }

      if(length(a) + length(b) > 0) {
        message(paste0("--", names(data)[i],": a total of ",length(a) + length(b), " features were removed due to the categories were not equal to 2!"))
      }
    }
    type[bindex] <- "binary"
  }

  fit <- LRAcluster_MOVICS(data, dimension = N.clust, types = as.list(type))
  dist <- fit$coordinate |> t() |> dist()
  clust.dend <- hclust(dist, method = clusterAlg)

  clustres <- data.frame(samID = colnames(data[[1]]),
                         clust = cutree(clust.dend,k = N.clust),
                         row.names = colnames(data[[1]]),
                         stringsAsFactors = FALSE)

  return(list(fit = fit, clust.res = clustres, clust.dend = clust.dend, mo.method = "LRAcluster"))
}

#MOVICS LRAcluster base code
#' Integrated analysis of cancer omics data by low-rank approximation using
#' LRAcluster code from the MOVICS package
#' @description
#' LRAcluster function from the MOVICS package. Used in the wrapper function
#' `getLRAcluster_MOVICS`.
#' @param data List of omics matrices, with samples in columns.
#' @param types List of data types, must be one of `"gaussian"` (default),
#' `"binomial"`, or `"poisson"`.
#' @param dimension The reduced dimension. Default is 2.
#' @param names Names of the datasets.
#' @returns A list of results. $coordinate, a matrix of the coordinates of the
#' samples in the reduced space. $potential, the ratio of explained variance.
#' @keywords internal
#' @references   * Lu, X., Meng, J., Zhou, Y., Jiang, L., and Yan, F. (2020).
#'   MOVICS: an R package for multi-omics integration and visualization in cancer
#'   subtyping. Bioinformatics, btaa1018.
#'   * Wu D, Wang D, Zhang MQ, Gu J (2015). Fast dimension reduction and
#'   integrative clustering of multi-omics data using low-rank approximation:
#'   application to cancer molecular classification. BMC Genomics, 16(1):1022.
LRAcluster_MOVICS <- function(data,
                              types,
                              dimension = 2,
                              names = as.character(1:length(data)))
{
  #--------#
  # binary #
  #--------#
  epsilon.binary<-2.0
  check.binary.row<-function(arr)
  {
    if (sum(!is.na(arr))==0)
    {
      return (F)
    }
    else
    {
      idx<-!is.na(arr)
      if (sum(arr[idx])==0 || sum(arr[idx])==sum(idx))
      {
        return (F)
      }
      else
      {
        return (T)
      }
    }
  }
  check.binary<-function(mat,name)
  {
    index<-apply(mat,1,check.binary.row)
    n<-sum(!index)
    if (n>0)
    {
      w<-paste("Warning: ",name," have ",as.character(n)," invalid lines!",sep="")
      warning(w)
    }
    mat_c<-mat[index,]
    rownames(mat_c)<-rownames(mat)[index]
    colnames(mat_c)<-colnames(mat)
    return (mat_c)
  }

  base.binary.row<-function(arr)
  {
    idx<-!is.na(arr)
    n<-sum(idx)
    m<-sum(arr[idx])
    return (log(m/(n-m)))
  }

  base.binary<-function(mat)
  {
    mat_b<-matrix(0,nrow(mat),ncol(mat))
    ar_b<-apply(mat,1,base.binary.row)
    mat_b[1:nrow(mat_b),]<-ar_b
    return (mat_b)
  }

  update.binary<-function(mat,mat_b,mat_now,eps)
  {
    mat_p<-mat_b+mat_now
    mat_u<-matrix(0,nrow(mat),ncol(mat))
    idx1<-!is.na(mat) & mat==1
    idx0<-!is.na(mat) & mat==0
    index<-is.na(mat)
    arr<-exp(mat_p)
    mat_u[index]<-mat_now[index]
    mat_u[idx1]<-mat_now[idx1]+eps*epsilon.binary/(1.0+arr[idx1])
    mat_u[idx0]<-mat_now[idx0]-eps*epsilon.binary*arr[idx0]/(1.0+arr[idx0])
    return (mat_u)
  }

  stop.binary<-function(mat,mat_b,mat_now,mat_u)
  {
    index<-!is.na(mat)
    mn<-mat_b+mat_now
    mu<-mat_b+mat_u
    arn<-exp(mn)
    aru<-exp(mu)
    idx1<-!is.na(mat) & mat==1
    idx0<-!is.na(mat) & mat==0
    lgn<-sum(log(arn[idx1]/(1+arn[idx1])))+sum(log(1/(1+arn[idx0])))
    lgu<-sum(log(aru[idx1]/(1+aru[idx1])))+sum(log(1/(1+aru[idx0])))
    return (lgu-lgn)
  }

  LL.binary<-function(mat,mat_b,mat_u)
  {
    index<-!is.na(mat)
    mu<-mat_b+mat_u
    aru<-exp(mu)
    idx1<-!is.na(mat) & mat==1
    idx0<-!is.na(mat) & mat==0
    lgu<-sum(log(aru[idx1]/(1+aru[idx1])))+sum(log(1/(1+aru[idx0])))
    return (lgu)
  }

  LLmax.binary<-function(mat)
  {
    return (0)
  }

  LLmin.binary<-function(mat,mat_b)
  {
    index<-!is.na(mat)
    aru<-exp(mat_b)
    idx1<-!is.na(mat) & mat==1
    idx0<-!is.na(mat) & mat==0
    lgu<-sum(log(aru[idx1]/(1+aru[idx1])))+sum(log(1/(1+aru[idx0])))
    return (lgu)
  }

  binary_type_base <- function( data,dimension=2 ,name="test")
  {
    data<-check.binary(data,name)
    data_b<-base.binary(data)
    data_now<-matrix(0,nrow(data),ncol(data))
    data_u<-update.binary(data,data_b,data_now)
    data_u<-nuclear_approximation(data_u,dimension)
    while (T)
    {
      thr<-stop.binary(data,data_b,data_now,data_u)
      message(thr)
      if (thr<0.2)
      {
        break
      }
      data_now<-data_u
      data_u<-update.binary(data,data_b,data_now)
      data_u<-nuclear_approximation(data_u,dimension)
    }
    return (data_now)
  }

  #----------#
  # gaussian #
  #----------#

  epsilon.gaussian=0.5

  check.gaussian.row<-function(arr)
  {
    if (sum(!is.na(arr))==0)
    {
      return (F)
    }
    else
    {
      return (T)
    }
  }
  check.gaussian<-function(mat,name)
  {
    index<-array(T,nrow(mat))
    for(i in 1:nrow(mat))
    {
      if (sum(is.na(mat[i,])==ncol(mat)))
      {
        war<-paste("Warning: ",name,"'s ",as.character(i)," line is all NA. Delete this line",sep="")
        warning(war)
        index[i]<-F
      }
    }
    mat_c<-mat[index,]
    rownames(mat_c)<-rownames(mat)[index]
    colnames(mat_c)<-colnames(mat)
    return (mat_c)
  }

  base.gaussian.row<-function(arr)
  {
    idx<-!is.na(arr)
    return (mean(arr[idx]))
  }

  base.gaussian<-function(mat)
  {
    mat_b<-matrix(0,nrow(mat),ncol(mat))
    ar_b<-apply(mat,1,base.gaussian.row)
    mat_b[1:nrow(mat_b),]<-ar_b
    return (mat_b)
  }

  update.gaussian<-function(mat,mat_b,mat_now,eps)
  {
    mat_p<-mat_b+mat_now
    mat_u<-matrix(0,nrow(mat),ncol(mat))
    index<-!is.na(mat)
    mat_u[index]<-mat_now[index]+eps*epsilon.gaussian*(mat[index]-mat_p[index])
    index<-is.na(mat)
    mat_u[index]<-mat_now[index]
    return (mat_u)
  }

  stop.gaussian<-function(mat,mat_b,mat_now,mat_u)
  {
    index<-!is.na(mat)
    mn<-mat_b+mat_now
    mu<-mat_b+mat_u
    ren<-mat[index]-mn[index]
    reu<-mat[index]-mu[index]
    lgn<- -0.5*sum(ren*ren)
    lgu<- -0.5*sum(reu*reu)
    return (lgu-lgn)
  }

  LL.gaussian<-function(mat,mat_b,mat_u)
  {
    index<-!is.na(mat)
    mu<-mat_b+mat_u
    reu<-mat[index]-mu[index]
    lgu<- -0.5*sum(reu*reu)
    return (lgu)
  }

  LLmax.gaussian<-function(mat)
  {
    return (0.0)
  }

  LLmin.gaussian<-function(mat,mat_b)
  {
    index<-!is.na(mat)
    reu<-mat[index]-mat_b[index]
    lgu<- -0.5*sum(reu*reu)
    return (lgu)
  }

  gaussian_base<-function(data,dimension=2,name="test")
  {
    data<-check.gaussian(data,name)
    data_b<-base.gaussian(data)
    data_now<-matrix(0,nrow(data),ncol(data))
    data_u<-update.gaussian(data,data_b,data_now)
    data_u<-nuclear_approximation(data_u,dimension)
    while(T)
    {
      thr<-stop.gaussian(data,data_b,data_now,data_u)
      message(thr)
      if (thr<0.2)
      {
        break
      }
      data_now<-data_u
      data_u<-update.gaussian(data,data_b,data_now)
      data_u<-nuclear_approximation(data_u,dimension)
    }
    return (data_now)
  }

  #---------#
  # poisson #
  #---------#

  epsilon.poisson<-0.5

  check.poisson.row<-function(arr)
  {
    if (sum(!is.na(arr))==0)
    {
      return (F)
    }
    else
    {
      idx<-!is.na(arr)
      if (sum(arr[idx]<0)>0)
      {
        return (F)
      }
      else
      {
        return (T)
      }
    }
  }

  check.poisson<-function(mat,name)
  {
    w<-paste(name," is poisson type. Add 1 to all counts",sep="")
    message(w)
    index<-apply(mat,1,check.poisson.row)
    n<-sum(!index)
    if (n>0)
    {
      w<-paste("Warning: ",name," have ",as.character(n)," invalid lines!",sep="")
      warning(w)
    }
    mat_c<-mat[index,]+1
    rownames(mat_c)<-rownames(mat)[index]
    colnames(mat_c)<-colnames(mat)
    return (mat_c)
  }

  base.poisson.row<-function(arr)
  {
    idx<-!is.na(arr)
    m<-sum(log(arr[idx]))
    n<-sum(idx)
    return(m/n)
  }

  base.poisson<-function(mat)
  {
    mat_b<-matrix(0,nrow(mat),ncol(mat))
    ar_b<-apply(mat,1,base.poisson.row)
    mat_b[1:nrow(mat_b),]<-ar_b
    return (mat_b)
  }

  update.poisson<-function(mat,mat_b,mat_now,eps)
  {
    mat_p<-mat_b+mat_now
    mat_u<-matrix(0,nrow(mat),ncol(mat))
    index<-!is.na(mat)
    mat_u[index]<-mat_now[index]+eps*epsilon.poisson*(log(mat[index])-mat_p[index])
    index<-is.na(mat)
    mat_u[index]<-mat_now[index]
    return (mat_u)
  }

  stop.poisson<-function(mat,mat_b,mat_now,mat_u)
  {
    index<-!is.na(mat)
    mn<-mat_b+mat_now
    mu<-mat_b+mat_u
    lgn<-sum(mat[index]*mn[index]-exp(mn[index]))
    lgu<-sum(mat[index]*mu[index]-exp(mu[index]))
    return (lgu-lgn)
  }

  LL.poisson<-function(mat,mat_b,mat_u)
  {
    index<-!is.na(mat)
    mu<-mat_b+mat_u
    lgu<-sum(mat[index]*mu[index]-exp(mu[index]))
    return (lgu)
  }

  LLmax.poisson<-function(mat)
  {
    index<-!is.na(mat)
    lgu<-sum(mat[index]*log(mat[index])-mat[index])
    return (lgu)
  }

  LLmin.poisson<-function(mat,mat_b)
  {
    index<-!is.na(mat)
    lgu<-sum(mat[index]*mat_b[index]-exp(mat_b[index]))
    return (lgu)
  }

  poisson_type_base<-function(data,dimension=2,name="test")
  {
    data<-check.poisson(data,name)
    data_b<-base.poisson(data)
    data_now<-matrix(0,nrow(data),ncol(data))
    data_u<-update.poisson(data,data_b,data_now)
    data_u<-nuclear_approximation(data_u,dimension)
    while(T)
    {
      thr<-stop.poisson(data,data_b,data_now,data_u)
      message(thr)
      if (thr<0.2)
      {
        break
      }
      data_now<-data_u
      data_u<-update.poisson(data,data_b,data_now)
      data_u<-nuclear_approximation(data_u,dimension)
    }
    return (data_now)
  }

  #----#
  # na #
  #----#

  nuclear_approximation<-function(mat,dimension)
  {
    svd<-svd(mat,nu=0,nv=0)
    if (dimension<length(svd$d))
    {
      lambda<-svd$d[dimension+1]
      svd<-svd(mat,nu=dimension,nv=dimension)
      indexh<-svd$d>lambda
      indexm<-svd$d<lambda
      dia<-array(svd$d,length(svd$d))
      dia[indexh]<-dia[indexh]-lambda
      dia[indexm]<-0
      mat_low<-svd$u%*%diag(c(dia[1:dimension],0))[1:dimension,1:dimension]%*%t(svd$v)
    }
    else
    {
      mat_low<-mat
    }
    return (mat_low)
  }

  #------------#
  # LRAcluster #
  #------------#
  check.matrix.element<-function(x)
  {
    if (!is.matrix(x))
    {
      return (T)
    }
    else
    {
      return (F)
    }
  }

  ncol.element<-function(x)
  {
    return (ncol(x))
  }

  nrow.element<-function(x)
  {
    return (nrow(x))
  }

  check<-function(mat,type,name)
  {
    if (type=="binary")
    {
      return (check.binary(mat,name))
    }
    else if (type=="gaussian")
    {
      return (check.gaussian(mat,name))
    }
    else if (type=="poisson")
    {
      return (check.poisson(mat,name))
    }
    else
    {
      e<-paste("unknown type ",type,sep="")
      stop(e)
    }
  }

  eps<-0.0
  if (!is.list(data))
  {
    stop("the input data must be a list!")
  }
  c<-sapply(data,check.matrix.element)
  if (sum(c)>0)
  {
    stop("each element of input list must be a matrix!")
  }
  c<-sapply(data,ncol.element)
  if (length(levels(factor(c)))>1)
  {
    stop("each element of input list must have the same column number!")
  }
  if (length(data)!=length(types))
  {
    stop("data and types must be the same length!")
  }
  nSample<-c[1]
  loglmin<-0
  loglmax<-0
  loglu<-0.0
  nData<-length(data)
  for (i in 1:nData)
  {
    data[[i]]<-check(data[[i]],types[[i]],names[[i]])
  }
  nGeneArr<-sapply(data,nrow.element)
  nGene<-sum(nGeneArr)
  indexData<-list()
  k=1
  for(i in 1:nData)
  {
    indexData[[i]]<- (k):(k+nGeneArr[i]-1)
    k<-k+nGeneArr[i]
  }
  base<-matrix(0,nGene,nSample)
  now<-matrix(0,nGene,nSample)
  update<-matrix(0,nGene,nSample)
  thr<-array(0,nData)
  for (i in 1:nData)
  {
    if (types[[i]]=="binary")
    {
      base[indexData[[i]],]<-base.binary(data[[i]])
      loglmin<-loglmin+LLmin.binary(data[[i]],base[indexData[[i]],])
      loglmax<-loglmax+LLmax.binary(data[[i]])
    }
    else if (types[[i]]=="gaussian")
    {
      base[indexData[[i]],]<-base.gaussian(data[[i]])
      loglmin<-loglmin+LLmin.gaussian(data[[i]],base[indexData[[i]],])
      loglmax<-loglmax+LLmax.gaussian(data[[i]])
    }
    else if (types[[i]]=="poisson")
    {
      base[indexData[[i]],]<-base.poisson(data[[i]])
      loglmin<-loglmin+LLmin.poisson(data[[i]],base[indexData[[i]],])
      loglmax<-loglmax+LLmax.poisson(data[[i]])
    }
  }
  for (i in 1:nData)
  {
    if (types[[i]]=="binary")
    {
      update[indexData[[i]],]<-update.binary(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
    }
    else if (types[[i]]=="gaussian")
    {
      update[indexData[[i]],]<-update.gaussian(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
    }
    else if (types[[i]]=="poisson")
    {
      update[indexData[[i]],]<-update.poisson(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
    }
  }
  update<-nuclear_approximation(update,dimension)
  nIter<-0
  thres<-array(Inf,3)
  epsN<-array(Inf,2)
  while(T)
  {
    for (i in 1:nData)
    {
      if (types[[i]]=="binary")
      {
        thr[i]<-stop.binary(data[[i]],base[indexData[[i]],],now[indexData[[i]],],update[indexData[[i]],])
      }
      else if (types[[i]]=="gaussian")
      {
        thr[i]<-stop.gaussian(data[[i]],base[indexData[[i]],],now[indexData[[i]],],update[indexData[[i]],])
      }
      else if (types[[i]]=="poisson")
      {
        thr[i]<-stop.poisson(data[[i]],base[indexData[[i]],],now[indexData[[i]],],update[indexData[[i]],])
      }
    }
    nIter<-nIter+1
    thres[1]<-thres[2]
    thres[2]<-thres[3]
    thres[3]<-sum(thr)
    epsN[1]<-epsN[2]
    epsN[2]<-eps
    if (nIter>5)
    {
      if (runif(1)<thres[1]*thres[3]/(thres[2]*thres[2]+thres[1]*thres[3]))
      {
        eps<-epsN[1]+0.05*runif(1)-0.025
      }
      else
      {
        eps<-epsN[2]+0.05*runif(1)-0.025
      }
      if (eps< -0.7)
      {
        eps<- 0
        epsN<-c(0,0)
      }
      if (eps > 1.4)
      {
        eps<-0
        epsN<-c(0,0)
      }
    }
    if (sum(thr)<nData*0.2)
    {
      break
    }
    now<-update
    for (i in 1:nData)
    {
      if (types[[i]]=="binary")
      {
        update[indexData[[i]],]<-update.binary(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
      }
      else if (types[[i]]=="gaussian")
      {
        update[indexData[[i]],]<-update.gaussian(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
      }
      else if (types[[i]]=="poisson")
      {
        update[indexData[[i]],]<-update.poisson(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
      }
    }
    update<-nuclear_approximation(update,dimension)
  }
  for (i in 1:nData)
  {
    if (types[[i]]=="binary")
    {
      loglu<-loglu+LL.binary(data[[i]],base[indexData[[i]],],update[indexData[[i]],])
    }
    else if (types[[i]]=="gaussian")
    {
      loglu<-loglu+LL.gaussian(data[[i]],base[indexData[[i]],],update[indexData[[i]],])
    }
    else if (types[[i]]=="poisson")
    {
      loglu<-loglu+LL.poisson(data[[i]],base[indexData[[i]],],update[indexData[[i]],])
    }
  }
  sv<-svd(update,nu=0,nv=dimension)
  coordinate<-diag(c(sv$d[1:dimension],0))[1:dimension,1:dimension]%*%t(sv$v)
  colnames(coordinate)<-colnames(data[[1]])
  rownames(coordinate)<-paste("PC ",as.character(1:dimension),sep="")
  ratio<-(loglu-loglmin)/(loglmax-loglmin)
  return (list("coordinate"=coordinate,"potential"=ratio))
}

#SNF affinity matrix creation per distance
#' Calculate affinity matrix as needed for SNF.
#' @description
#' Wrapper function to calculate the affinity matrix needed for SNF. All
#' distance metrics provided in the EINS method can be used.
#' @param Distance String indicating the distance used to calculate the
#' distance matrix. Must be one of: `"euclidean squared"`, `"euclidean"`,
#' `"manhattan"`, `"minkowski 0.25"`, `"minkowski 0.5"`, `"minkowski 3"`,
#' `"minkowski 4"`.
#' @param Data List of omics matrices, with samples in rows.
#' @param Neighbours Number of nearest neighbors.
#' @param Sigma Variance for the local model.
#' @returns Affinity matrix calculated with the distance metric provided.
#' @keywords internal
#' @references Wang, B., Mezlini, A., Demir, F., Fiume, M., Tu, Z., Brudno,
#' M., Haibe-Kains, B., and Goldenberg, A. (2021) SNFtool: Similarity
#' Network Fusion. R package version 2.3.1
SNF_Distance_Affinity_Matrix = function(Distance = NULL,
                                        Data = NULL,
                                        Neighbours = NULL,
                                        Sigma = NULL){
  if(Distance == "euclidean squared"){
    AffMatrix = lapply(Data, function(dat){
      dist = SNFtool::dist2(dat, dat)
      sim = SNFtool::affinityMatrix(dist,
                                    K = Neighbours,
                                    sigma = Sigma)
    })
  } else if(Distance == "euclidean"){
    AffMatrix = lapply(Data, function(dat){
      dist = as.matrix(stats::dist(dat, method = "euclidean"))
      sim = SNFtool::affinityMatrix(dist,
                                    K = Neighbours,
                                    sigma = Sigma)
    })
  } else if(Distance == "manhattan"){
    AffMatrix = lapply(Data, function(dat){
      dist = as.matrix(stats::dist(dat, method = "manhattan"))
      sim = SNFtool::affinityMatrix(dist,
                                    K = Neighbours,
                                    sigma = Sigma)
    })
  } else if(Distance == "minkowski 0.25"){
    AffMatrix = lapply(Data, function(dat){
      dist = as.matrix(stats::dist(dat, method = "minkowski", p = 0.25))
      sim = SNFtool::affinityMatrix(dist,
                                    K = Neighbours,
                                    sigma = Sigma)
    })
  } else if(Distance == "minkowski 0.5"){
    AffMatrix = lapply(Data, function(dat){
      dist = as.matrix(stats::dist(dat, method = "minkowski", p = 0.5))
      sim = SNFtool::affinityMatrix(dist,
                                    K = Neighbours,
                                    sigma = Sigma)
    })
  } else if(Distance == "minkowski 3"){
    AffMatrix = lapply(Data, function(dat){
      dist = as.matrix(stats::dist(dat, method = "minkowski", p = 3))
      sim = SNFtool::affinityMatrix(dist,
                                    K = Neighbours,
                                    sigma = Sigma)
    })
  } else if(Distance == "minkowski 4"){
    AffMatrix = lapply(Data, function(dat){
      dist = as.matrix(stats::dist(dat, method = "minkowski", p = 4))
      sim = SNFtool::affinityMatrix(dist,
                                    K = Neighbours,
                                    sigma = Sigma)
    })
  } else {stop("Unknown distance measure selected. Choose one of the following options: euclidean squared (default), euclidean, manhattan, minkowski 0.25, minkowski 0.5, minkowski 3, minkowski 4")}

  return(AffMatrix)
}

#Min Max scaling
#' Min max scaling
#' @description
#' Min max scaling into scale 0-1. Used in function `Feature_Weight_Scaling`.
#' @param x Dataframe with single numerical column.
#' @returns Scaled dataframe.
#' @keywords internal
OverlapScaling = function(x){
  if(min(x) < 0){
    2 * ((x - min(x)) / (max(x) - min(x))) - 1
  }else{
    (x - min(x)) / (max(x) - min(x))
  }
}

#scale feature data for factors
#' Scale feature weight data from multi-omics integration
#' @description
#' Scale the feature weight data calculated by the multi-omics integration
#' methods, per method, number of factors/components and omics type. Scale used
#' is 0-1. Used in EINS R6 class functions `run_Multi_Omics_Feature_HClust` and
#' `run_Feature_Ensemble_HClust`.
#' @param FeatureData List of feature weight dataframes for the different
#' multi-omics integration methods, by number of factors/components and omics
#' type.
#' @returns List of scaled feature weight dataframes, stored by multi-omics
#' integration method, number of factors/components and omics type.
#' @keywords internal
Feature_Weight_Scaling = function(FeatureData){
  ScaledFeatures = list()
  FeatMethods = names(FeatureData)
  for(method in FeatMethods){
    FeatFactors = names(FeatureData[[method]])
    for(factor in FeatFactors){
      FeatOmics = names(FeatureData[[method]][[factor]])
      for(omics in FeatOmics){
        Data = FeatureData[[method]][[factor]][[omics]]
        Data$Features = as.character(Data$Features)
        Data = Data |>
          dplyr::rowwise() |>
          dplyr::mutate(Features = dplyr::case_when(
            (strsplit(Features, split = "_")[[1]][2]) == omics ~ Features,
            TRUE ~ paste0(Features, "_", omics)))
        Data = tibble::column_to_rownames(Data, var = "Features")
        Data[sapply(Data, is.character)] = lapply(Data[sapply(Data, is.character)], as.numeric)
        ScaleData = OverlapScaling(Data)
        ScaleData = tibble::rownames_to_column(ScaleData, var = "Features")
        ScaledFeatures[[method]][[factor]][[omics]] = ScaleData
      }
    }
  }
  return(ScaledFeatures)
}

#feature clustering per multi-omics method
#' Cluster scaled feature weight data from multi-omics integration
#' @description
#' Use scaled feature weight data per multi-omics integration method as
#' calculated by `Feature_Weight_Scaling()` to calculate a feature distance
#' matrix per multi-omics integration method. Hierarchical clustering is then
#' performed on these matrices to create a dendrogram of features per
#' multi-omics integration method.
#' @param ScaledFeaturesData List of scaled feature weight dataframes for the
#' different multi-omics integration methods, by number of factors/components
#' and omics type, as calculated with `Feature_Weight_Scaling()`.
#' @param nFactors Number of factors (components) for which the hierarchical
#' clustering is to be performed. All multi-omics integration methods will need
#' to have been performed with this number of factors/components.
#' @param Distance Distance metric to be used for the calculation of the
#' feature distance matrices. Must be one of: `"euclidean"` (default),
#' `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"` or `"minkowski"`.
#' If `"minkowski"` is selected, argument `"MinkowskiPower"` needs to be
#' included.
#' @param MinkowskiPower Power of the Minkowski distance. Default is `NULL`.
#' @param Linkage Agglomeration method to be used for the hierarchical
#' clustering. Must be one of: `"ward.D"`, `"ward.D2"` (default), `"single"`,
#' `"complete"`, `"average"`, `"mcquitty"`, `"median"` or `"centroid"`.
#' @returns A list of results. Feature distance matrices per multi-omics
#' integration method, stored in $DistMat. Hierarchical clustering trees per
#' multi-omics integration method, stored in $HClust. Feature dendrogram per
#' multi-omics integration method, stored in $Dendrogram.
#' @keywords internal
Multi_Omics_Feature_Clustering = function(ScaledFeaturesData,
                                          nFactors = NULL,
                                          Distance = "euclidean",
                                          MinkowskiPower = NULL,
                                          Linkage = "ward.D2"){
  NamedFact = list()
  MethodComb = list()
  OmicsComb = list()
  MethodFeatWeight = list()
  FeatMethods = names(ScaledFeaturesData)
  for(method in FeatMethods){
    FeatFactors = names(ScaledFeaturesData[[method]])
    for(factor in FeatFactors){
      n = as.numeric(strsplit(factor, split = "_")[[1]][2])
      if(n == nFactors){
        FeatOmics = names(ScaledFeaturesData[[method]][[factor]])
        for(omics in FeatOmics){
          Data = subset(ScaledFeaturesData[[method]][[factor]][[omics]], select = 1:(n+1))
          for(i in 1:n){
            Name = paste0(method, "_F", i)
            colnames(Data)[i+1] = Name
          }
          rownames(Data) = NULL
          NamedFact[[method]][[factor]][[omics]] = Data
          MethodComb[[method]][[factor]] = data.table::rbindlist(NamedFact[[method]][[factor]], fill = TRUE)
          OmicsComb[[method]] = purrr::reduce(MethodComb[[method]], dplyr::full_join, by = "Features")
          MethodFeatWeight[[method]] = tibble::column_to_rownames(OmicsComb[[method]], var = "Features")
          MethodFeatWeight[[method]] = t(as.matrix(MethodFeatWeight[[method]]))
        }
      }
    }
  }

  DistList = list()
  for(method in FeatMethods){
    Data = MethodFeatWeight[[method]]
    tData = t(Data)
    #dist creates distance matrix for rows
    if(Distance == "minkowski"){
      if(is.null(MinkowskiPower)){
        DistList[[method]] = stats::dist(tData,
                                         method = Distance)
      }else{
        DistList[[method]] = stats::dist(tData,
                                         method = Distance,
                                         p = MinkowskiPower)
      }
    }else{
      DistList[[method]] = stats::dist(tData,
                                       method = Distance)
    }
  }

  HClustList = list()
  for(method in FeatMethods){
    HClustList[[method]] = fastcluster::hclust(d = DistList[[method]], method = Linkage)
  }

  DendroList = list()
  for(method in FeatMethods){
    DendroList[[method]] = stats::as.dendrogram(HClustList[[method]])
  }

  MOHClust = list()
  MOHClust$DistMat = DistList
  MOHClust$HClust = HClustList
  MOHClust$Dendrogram = DendroList
  return(MOHClust)
}

#----------------Ensemble integration----------------
#Create list of cluster results for ensemble clustering
#' Prepare the cluster results from the different multi-omics integration
#' methods for ensemble clustering
#' @description
#' Create a list of dataframes with the sample clustering results of the
#' multi-omics integration methods with all calculated cluster numbers. This
#' dataframe is used for ensemble clustering in `Ensemble_Cluster_CHC()`.
#' @param Data A list of dataframes, stored by multi-omics integration method
#' and cluster number. Each dataframe contains a single cluster result.
#' @returns A list of dataframes with samples in rows and cluster assignment
#' in columns per multi-omics integration method.
#' @keywords internal
Ensemble_Cluster_Data = function(Data){
  AllClusterRes = list()

  ClusterMethods = names(Data)
  for(method in ClusterMethods){
    if(method == "SNF"){
      ClusterSNFDistance = names(Data$SNF)
      for(distance in ClusterSNFDistance){
        ClusterSNFClustnum = names(Data$SNF[[distance]])
        for(clustnum in ClusterSNFClustnum){
          SNFClusters = as.integer(strsplit(clustnum, "_")[[1]][2])
          SNFName = paste0(method, "_", distance, "_", SNFClusters)
          AllClusterRes[[SNFName]] = Data$SNF[[distance]][[clustnum]]
        }
      }
    }else if(method == "GAUDI"){
      ClusterGAUDIClustmethod = names(Data$GAUDI)
      for(clustmethod in ClusterGAUDIClustmethod){
        ClusterGAUDIClustnum = names(Data$GAUDI[[clustmethod]])
        for(clustnum in ClusterGAUDIClustnum){
          GAUDIClusters = as.integer(strsplit(clustnum, "_")[[1]][2])
          GAUDIName = paste0(method, "_", clustmethod, "_", GAUDIClusters)
          AllClusterRes[[GAUDIName]] = Data$GAUDI[[clustmethod]][[clustnum]]
        }
      }
    }else{
      ClusterClustnum = names(Data[[method]])
      for(clustnum in ClusterClustnum){
        MethodClusters = as.integer(strsplit(clustnum, "_")[[1]][2])
        MethodName = paste0(method, "_", MethodClusters)
        AllClusterRes[[MethodName]] = Data[[method]][[clustnum]]
      }
    }
  }
  return(AllClusterRes)
}

#Sample co-association matrix from a set of cluster assignments
#' Build the sample co-association matrix from a set of cluster assignments
#' @description
#' For every pair of samples, the fraction of the supplied cluster assignments
#' that place them in the same cluster. This is the consensus matrix that
#' `Ensemble_Cluster_CHC()`, `Ensemble_Cluster_MDS()`,
#' `Ensemble_Cluster_MDS_HC()`and `Ensemble_Cluster_MPA()` all operate on,
#' factored out so that every method consumes an identical matrix rather than
#' rebuilding it.
#' @param MethodClusterResults A list of dataframes with cluster assignments per
#' multi-omics integration method and cluster number, as created by
#' `Ensemble_Cluster_Data()`.
#' @param Sample_IDs Vector of sample names.
#' @returns A square matrix of co-association frequencies averaged across
#' multi-omics integration sample clustering results, samples in rows and columns,
#' values between 0 and 1.
#' @keywords internal
#' @references Lu X, Meng J, Zhou Y, Jiang L, Yan F. MOVICS: an R package for
#' multi-omics integration and visualization in cancer subtyping. Bioinformatics.
#' 2021 Apr 1;36(22-23):5539-5541.
Ensemble_CoAssociation_Matrix = function(MethodClusterResults,
                                         Sample_IDs){
  Ens_Matrix = matrix(0, nrow = length(Sample_IDs), ncol = length(Sample_IDs))
  for(i in 1:length(MethodClusterResults)){
    Ens_Clust_Result = MethodClusterResults[[i]]
    Ens_Clust_Group = Ens_Clust_Result$Cluster
    names(Ens_Clust_Group) = Sample_IDs
    Ens_Clust_Ans = matrix(0, nrow = length(Sample_IDs), ncol = length(Sample_IDs))
    rownames(Ens_Clust_Ans) = colnames(Ens_Clust_Ans) = Sample_IDs
    Ens_Clust_Unique_Group = unique(Ens_Clust_Group)

    for(j in 1:length(Ens_Clust_Unique_Group)){
      Ens_Clust_Group_j = names(Ens_Clust_Group)[Ens_Clust_Group == Ens_Clust_Unique_Group[j]]
      Ens_Clust_Ans[Ens_Clust_Group_j, Ens_Clust_Group_j] = 1
    }
    Ens_Matrix = Ens_Matrix + as.matrix(Ens_Clust_Ans)
  }
  Ens_Sim_Matrix = Ens_Matrix/length(MethodClusterResults)
  rownames(Ens_Sim_Matrix) = colnames(Ens_Sim_Matrix) = Sample_IDs
  return(Ens_Sim_Matrix)
}

#Create list of cluster results for co-association clustering
#' Prepare the cluster results from the different multi-omics integration
#' methods for co-association clustering
#' @description
#' Create a list of dataframes with the sample clustering results of the
#' multi-omics integration methods with the specified cluster number. This
#' dataframe is used for ensemble clustering in `Ensemble_Cluster_CHC()`.
#' @param Data A list of dataframes, stored by multi-omics integration method
#' and cluster number. Each dataframe contains a single cluster result.
#' @returns A list of dataframes with samples in rows and cluster assignment
#' in columns per multi-omics integration method.
#' @keywords internal
Ensemble_Cluster_Data_CCC = function(Data,
                                    ResClust = NULL){
  AllClusterRes = list()

  ClusterMethods = names(Data)
  for(method in ClusterMethods){
    if(method == "SNF"){
      ClusterSNFDistance = names(Data$SNF)
      for(distance in ClusterSNFDistance){
        ClusterSNFClustnum = names(Data$SNF[[distance]])
        for(clustnum in ClusterSNFClustnum){
          SNFClusters = as.integer(strsplit(clustnum, "_")[[1]][2])
          if(SNFClusters == ResClust){
            SNFName = paste0(method, "_", distance, "_", SNFClusters)
            AllClusterRes[[SNFName]] = Data$SNF[[distance]][[clustnum]]
          }
        }
      }
    }else if(method == "GAUDI"){
      ClusterGAUDIClustmethod = names(Data$GAUDI)
      for(clustmethod in ClusterGAUDIClustmethod){
        ClusterGAUDIClustnum = names(Data$GAUDI[[clustmethod]])
        for(clustnum in ClusterGAUDIClustnum){
          GAUDIClusters = as.integer(strsplit(clustnum, "_")[[1]][2])
          if(GAUDIClusters == ResClust){
            GAUDIName = paste0(method, "_", clustmethod, "_", GAUDIClusters)
            AllClusterRes[[GAUDIName]] = Data$GAUDI[[clustmethod]][[clustnum]]
          }
        }
      }
    }else{
      ClusterClustnum = names(Data[[method]])
      for(clustnum in ClusterClustnum){
        MethodClusters = as.integer(strsplit(clustnum, "_")[[1]][2])
        if(MethodClusters == ResClust){
          MethodName = paste0(method, "_", MethodClusters)
          AllClusterRes[[MethodName]] = Data[[method]][[clustnum]]
        }
      }
    }
  }
  return(AllClusterRes)
}

#Create ensemble clustering from ensemble matrix
#' Create an ensemble clustering assignment from all available multi-omics
#' integration results
#' @description
#' Using the dataframe with multi-omics integration cluster assignments as
#' created in `Ensemble_Cluster_Data()`, create a sample similarity matrix by
#' counting the number of times samples cluster together. Code for similarity
#' matrix is adapted from the MOVICS pacakage. Hierarchical
#' clustering is then performed on this matrix to provide a more robust
#' clustering result by combining multiple cluster assignments.
#' @param MethodClusterResults A list of dataframes with cluster assignments per
#' multi-omics integration method.
#' @param Sample_IDs Vector of sample names.
#' @param Distance Distance metric to be used for the calculation of the
#' feature distance matrices. Must be one of: `"euclidean"` (default),
#' `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"` or `"minkowski"`.
#' If `"minkowski"` is selected, argument `"MinkowskiPower"` needs to be
#' included.
#' @param MinkowskiPower Power of the Minkowski distance. Default is `NULL`.
#' @param Linkage Agglomeration method to be used for the hierarchical
#' clustering. Must be one of: `"ward.D"`, `"ward.D2"` (default), `"single"`,
#' `"complete"`, `"average"`, `"mcquitty"`, `"median"` or `"centroid"`.
#' @param Clusters Number of clusters to create from the hierarchical clustering
#' result. Default is `NULL`.
#' @returns A list of results. Cluster assignment from the hierarchical
#' clustering result, into the number of cluster as provided in `Clusters`,
#' stored in $ClusterRes. Sample distance matrix, stored in $DistMat.
#' hierarchical clustering tree, stored in $HClustRes.
#' @keywords internal
#' @references Lu X, Meng J, Zhou Y, Jiang L, Yan F. MOVICS: an R package for
#' multi-omics integration and visualization in cancer subtyping. Bioinformatics.
#' 2021 Apr 1;36(22-23):5539-5541.
Ensemble_Cluster_CHC = function(MethodClusterResults,
                                Sample_IDs,
                                Distance = "euclidean",
                                MinkowskiPower = NULL,
                                Linkage = "ward.D2",
                                Clusters = NULL){
  Ens_Sim_Matrix = Ensemble_CoAssociation_Matrix(MethodClusterResults = MethodClusterResults,
                                                 Sample_IDs = Sample_IDs)
  Ens_Dissim_Matrix = 1 - Ens_Sim_Matrix
  if(Distance == "minkowski"){
    if(is.null(MinkowskiPower)){
      Dist_Matrix = stats::dist(Ens_Dissim_Matrix,
                                method = Distance)
    }else{
      Dist_Matrix = stats::dist(Ens_Dissim_Matrix,
                                method = Distance,
                                p = MinkowskiPower)
    }
  }else{
    Dist_Matrix = stats::dist(Ens_Dissim_Matrix,
                              method = Distance)
  }
  Ens_HClust = fastcluster::hclust(Dist_Matrix, method = Linkage)
  Ens_Clust_Result = stats::cutree(Ens_HClust, k = Clusters)

  EnsembleClustering = list()
  EnsembleClustering$ClusterRes = Ens_Clust_Result
  EnsembleClustering$DistMat = Dist_Matrix
  EnsembleClustering$HClustRes = Ens_HClust

  return(EnsembleClustering)
}

#Create ensemble clustering from ensemble matrix
#' Create an ensemble clustering assignment from all available multi-omics
#' integration results and perform Multidimensional Scaling (MDS) followed
#' by hierarchical clustering
#' @description
#' Using the dataframe with multi-omics integration cluster assignments as
#' created in `Ensemble_Cluster_Data()`, create a sample similarity matrix by
#' counting the number of times samples cluster together. Code for similarity
#' matrix is adapted from the MOVICS package. Classical MDS is
#' applied to this matrix to obtain sample coordinates. Hierarchical
#' clustering is then performed on this MDS matrix. The MDS dimension is
#' selected automatically, by using leave-one-method-out stabilities, to favor
#' the most compact representation which preserves robustness. Stability is
#' quantified as the Adjusted Rand Index  between the full-ensemble MDS-HC
#' clustering and each leave-one-method-out MDS-HC clustering.
#' @param MethodClusterResults A list of dataframes with cluster assignments per
#' multi-omics integration method.
#' @param Sample_IDs Vector of sample names.
#' @param Distance Distance metric to be used for the calculation of the
#' feature distance matrices. Must be one of: `"euclidean"` (default),
#' `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"` or `"minkowski"`.
#' If `"minkowski"` is selected, argument `"MinkowskiPower"` needs to be
#' included.
#' @param MinkowskiPower Power of the Minkowski distance. Default is `NULL`.
#' @param CandidateDimensions Number of MDS dimensions to be tested. Default is
#' 2:10 dimensions.
#' @param StabilityEpsilon How much the selected leave-one-method-out run may
#' diverge from the maximum observed stability, default is 0.02.
#' @param AddConstant Logical, indicating if an additive constant c* should
#' be computed, and added to the non-diagonal dissimilarities such that the
#' modified dissimilarities are Euclidean.
#' @param Linkage Agglomeration method to be used for the hierarchical
#' clustering. Must be one of: `"ward.D"`, `"ward.D2"` (default), `"single"`,
#' `"complete"`, `"average"`, `"mcquitty"`, `"median"` or `"centroid"`.
#' @param Clusters Number of clusters to create from the hierarchical clustering
#' result. Default is `NULL`.
#' @returns A list of results. Sample clustering results, stored in $ClusterRes.
#' Sample distance matrix, stored in $DistMat. Consensus matrix, stored in
#' $ConsensusMatrix. Dissimilarity matrix, stored in $DissimilarityMatrix. MDS
#' final embedding, stored in $Embed. MDS embedding distance matrix, $EmbedDist.
#' Full final MDS result, provided in $MDS. Dimension calculation information,
#' provided in $DimensionInfo. Hierarchical clustering tree, stored in $HClustRes.
#' @keywords internal
#' @references Lu X, Meng J, Zhou Y, Jiang L, Yan F. MOVICS: an R package for
#' multi-omics integration and visualization in cancer subtyping. Bioinformatics.
#' 2021 Apr 1;36(22-23):5539-5541.
Ensemble_Cluster_MDS_HC = function(MethodClusterResults,
                                   Sample_IDs,
                                   Distance = "euclidean",
                                   MinkowskiPower = NULL,
                                   CandidateDimensions = c(2:10),
                                   StabilityEpsilon = 0.02,
                                   AddConstant = TRUE,
                                   Linkage = "ward.D2",
                                   Clusters = NULL){
  if(is.null(Clusters)){
    stop("Number of clusters must be provided.")
  }

  ## Internal helper

  #Adjusted Rand Index
  AdjustedRandIndex = function(x, y){
    x = as.factor(x)
    y = as.factor(y)

    tab = table(x, y)

    comb2 = function(z){
      z * (z - 1) / 2
    }

    sum_comb = sum(comb2(tab))
    sum_row = sum(comb2(rowSums(tab)))
    sum_col = sum(comb2(colSums(tab)))

    n = sum(tab)
    total_comb = comb2(n)

    expected_index = (sum_row * sum_col) / total_comb
    max_index = 0.5 * (sum_row + sum_col)

    if(max_index == expected_index){
      return(0)
    }

    ari = (sum_comb - expected_index) / (max_index - expected_index)
    return(as.numeric(ari))
  }

  #extract cluster vector
  GetClusterVector = function(Ens_Clust_Result,
                              Sample_IDs){
    if(is.data.frame(Ens_Clust_Result) || is.matrix(Ens_Clust_Result)){
      if("Cluster" %in% colnames(Ens_Clust_Result)){
        Ens_Clust_Group = Ens_Clust_Result[, "Cluster"]
      }else{
        Ens_Clust_Group = Ens_Clust_Result[, 1]
      }

      if(!is.null(rownames(Ens_Clust_Result)) && all(Sample_IDs %in% rownames(Ens_Clust_Result))){
        names(Ens_Clust_Group) = rownames(Ens_Clust_Result)
        Ens_Clust_Group = Ens_Clust_Group[Sample_IDs]
      }else{
        names(Ens_Clust_Group) = Sample_IDs
      }

    }else if(is.list(Ens_Clust_Result) && !is.null(Ens_Clust_Result$Cluster)){
      Ens_Clust_Group = Ens_Clust_Result$Cluster
      if(!is.null(names(Ens_Clust_Group)) && all(Sample_IDs %in% names(Ens_Clust_Group))){
        Ens_Clust_Group = Ens_Clust_Group[Sample_IDs]
      }else{
        names(Ens_Clust_Group) = Sample_IDs
      }
    }else{
      Ens_Clust_Group = Ens_Clust_Result
      if(!is.null(names(Ens_Clust_Group)) && all(Sample_IDs %in% names(Ens_Clust_Group))){
        Ens_Clust_Group = Ens_Clust_Group[Sample_IDs]
      }else{
        names(Ens_Clust_Group) = Sample_IDs
      }
    }
    return(Ens_Clust_Group)
  }

  # Build consensus matrix and dissimilarity matrix
  BuildConsensus = function(MethodClusterResults,
                            Sample_IDs,
                            DropIndex = NULL){
    if(!is.null(DropIndex)){
      MethodClusterResults = MethodClusterResults[-DropIndex]
    }

    if(length(MethodClusterResults) < 1){
      stop("Too few clustering results available to perform stability selection, perform additional multi-omics integration methods to calculate ensemble sample MDS_HC clustering result.")
    }

    Ens_Clust_Samples = Sample_IDs
    N_Samples = length(Ens_Clust_Samples)

    Ens_Matrix = matrix(0,
                        nrow = N_Samples,
                        ncol = N_Samples)

    rownames(Ens_Matrix) = colnames(Ens_Matrix) = Ens_Clust_Samples

    for(i in seq_along(MethodClusterResults)){
      Ens_Clust_Result = MethodClusterResults[[i]]
      Ens_Clust_Group = GetClusterVector(Ens_Clust_Result = Ens_Clust_Result,
                                         Sample_IDs = Ens_Clust_Samples)

      Ens_Clust_Ans = matrix(0,
                             nrow = N_Samples,
                             ncol = N_Samples)

      rownames(Ens_Clust_Ans) = colnames(Ens_Clust_Ans) = Ens_Clust_Samples
      Ens_Clust_Unique_Group = unique(Ens_Clust_Group)

      for(j in seq_along(Ens_Clust_Unique_Group)){
        Ens_Clust_Group_j = names(Ens_Clust_Group)[Ens_Clust_Group == Ens_Clust_Unique_Group[j]]
        Ens_Clust_Ans[Ens_Clust_Group_j, Ens_Clust_Group_j] = 1
      }
      Ens_Matrix = Ens_Matrix + as.matrix(Ens_Clust_Ans)
    }
    Ens_Sim_Matrix = Ens_Matrix / length(MethodClusterResults)
    rownames(Ens_Sim_Matrix) = colnames(Ens_Sim_Matrix) = Ens_Clust_Samples
    Ens_Dissim_Matrix = 1 - Ens_Sim_Matrix
    diag(Ens_Dissim_Matrix) = 0
    rownames(Ens_Dissim_Matrix) = colnames(Ens_Dissim_Matrix) = Ens_Clust_Samples
    Dist_Matrix = stats::as.dist(Ens_Dissim_Matrix)

    return(list(ConsensusMatrix = Ens_Sim_Matrix,
                DissimilarityMatrix = Ens_Dissim_Matrix,
                DistMat = Dist_Matrix))
  }
  #cluster MDS points using HC
  ClusterFromPoints = function(MDS_Points,
                               Sample_IDs,
                               Distance = "euclidean",
                               MinkowskiPower = NULL,
                               Linkage = "ward.D2",
                               Clusters = NULL){
    if(is.null(Clusters)){
      stop("Number of Clusters must be provided.")
    }

    MDS_Points = as.matrix(MDS_Points)
    rownames(MDS_Points) = Sample_IDs

    if(Distance == "minkowski"){
      if(is.null(MinkowskiPower)){
        Embed_Dist_Matrix = stats::dist(MDS_Points,
                                        method = Distance)
      }else{
        Embed_Dist_Matrix = stats::dist(MDS_Points,
                                        method = Distance,
                                        p = MinkowskiPower)
      }
    }else{
      Embed_Dist_Matrix = stats::dist(MDS_Points,
                                      method = Distance)
    }

    HClust_Result = stats::hclust(d = Embed_Dist_Matrix,
                                  method = Linkage)

    Cluster_Result = stats::cutree(tree = HClust_Result,
                                   k = Clusters)

    Cluster_Result = Cluster_Result[Sample_IDs]

    return(list(ClusterRes = Cluster_Result,
                EmbedDist = Embed_Dist_Matrix,
                HClustRes = HClust_Result))
  }

  ## Main function body
  Ens_Clust_Samples = Sample_IDs
  N_Samples = length(Ens_Clust_Samples)

  if(length(MethodClusterResults) < 3){
    stop("MDS_HC stability selection requires at least three multi-omics integration methods for sample clustering. Run additional multi-omics integration methods to perform `run_Ensemble_Sample_MDS_HC()`.")
  }

  Max_Possible_Dimensions = N_Samples - 1
  CandidateDimensions = sort(unique(as.integer(CandidateDimensions)))
  CandidateDimensions = CandidateDimensions[CandidateDimensions >= 1 & CandidateDimensions <= Max_Possible_Dimensions]

  if(length(CandidateDimensions) == 0){
    stop("No valid CandidateDimensions available. CandidateDimensions must fall between 1 and n Samples - 1.")
  }

  Max_Candidate_Dim = max(CandidateDimensions)

  ## Full ensemble consensus and full MDS
  Full_Consensus = BuildConsensus(MethodClusterResults = MethodClusterResults,
                                  Sample_IDs = Ens_Clust_Samples)

  Full_MDS = stats::cmdscale(d = Full_Consensus$DistMat,
                             k = Max_Candidate_Dim,
                             eig = TRUE,
                             add = AddConstant)

  Full_Points = as.matrix(Full_MDS$points)
  rownames(Full_Points) = Ens_Clust_Samples

  Full_Labels = list()

  for(p in CandidateDimensions){
    Full_Labels[[paste0("Dim_", p)]] = ClusterFromPoints(MDS_Points = Full_Points[, seq_len(p), drop = FALSE],
                                                         Sample_IDs = Ens_Clust_Samples,
                                                         Distance = Distance,
                                                         MinkowskiPower = MinkowskiPower,
                                                         Linkage = Linkage,
                                                         Clusters = Clusters)$ClusterRes
  }
  ## Leave-one-method-out stability selection
  N_Methods = length(MethodClusterResults)

  Stability_Matrix = matrix(NA,
                            nrow = length(CandidateDimensions),
                            ncol = N_Methods)

  rownames(Stability_Matrix) = paste0("Dim_", CandidateDimensions)
  colnames(Stability_Matrix) = paste0("Drop_", seq_len(N_Methods))

  for(drop_i in seq_len(N_Methods)){
    Drop_Consensus = BuildConsensus(MethodClusterResults = MethodClusterResults,
                                    Sample_IDs = Ens_Clust_Samples,
                                    DropIndex = drop_i)

    Drop_MDS = stats::cmdscale(d = Drop_Consensus$DistMat,
                               k = Max_Candidate_Dim,
                               eig = TRUE,
                               add = AddConstant)

    Drop_Points = as.matrix(Drop_MDS$points)
    rownames(Drop_Points) = Ens_Clust_Samples
    for(p in CandidateDimensions){
      key = paste0("Dim_", p)

      Drop_Labels = ClusterFromPoints(MDS_Points = Drop_Points[, seq_len(p), drop = FALSE],
                                      Sample_IDs = Ens_Clust_Samples,
                                      Distance = Distance,
                                      MinkowskiPower = MinkowskiPower,
                                      Linkage = Linkage,
                                      Clusters = Clusters)$ClusterRes

      Stability_Matrix[key, paste0("Drop_", drop_i)] = AdjustedRandIndex(Full_Labels[[key]],
                                                                         Drop_Labels)
    }
  }

  Mean_Stability = rowMeans(Stability_Matrix,
                            na.rm = TRUE)

  SD_Stability = apply(Stability_Matrix,
                       1,
                       stats::sd,
                       na.rm = TRUE)

  Max_Stability = max(Mean_Stability,
                      na.rm = TRUE)

  Eligible = CandidateDimensions[Mean_Stability >= (Max_Stability - StabilityEpsilon)]

  SelectedDimensions = min(Eligible)

  Stability_Table = data.frame(MDSDimensions = CandidateDimensions,
                               MeanARI = as.numeric(Mean_Stability),
                               SDARI = as.numeric(SD_Stability),
                               stringsAsFactors = FALSE)

  ## Final MDS-HC using selected dimensions on full consensus matrix
  Final_MDS = stats::cmdscale(d = Full_Consensus$DistMat,
                              k = SelectedDimensions,
                              eig = TRUE,
                              add = AddConstant)
  Final_MDS_Points = as.matrix(Final_MDS$points)
  rownames(Final_MDS_Points) = Ens_Clust_Samples
  colnames(Final_MDS_Points) = paste0("MDS", seq_len(ncol(Final_MDS_Points)))

  Final_Clustering = ClusterFromPoints(MDS_Points = Final_MDS_Points,
                                       Sample_IDs = Ens_Clust_Samples,
                                       Distance = Distance,
                                       MinkowskiPower = MinkowskiPower,
                                       Linkage = Linkage,
                                       Clusters = Clusters)

  DimensionInfo = list()
  DimensionInfo$SelectionMethod = "stability"
  DimensionInfo$SelectedDimensions = SelectedDimensions
  DimensionInfo$CandidateDimensions = CandidateDimensions
  DimensionInfo$StabilityEpsilon = StabilityEpsilon
  DimensionInfo$StabilityMatrix = Stability_Matrix
  DimensionInfo$StabilityTable = Stability_Table
  DimensionInfo$MaxStability = Max_Stability
  DimensionInfo$FullMDSForSelection = Full_MDS
  DimensionInfo$AddConstant = AddConstant
  DimensionInfo$Linkage = Linkage
  DimensionInfo$Distance = Distance
  DimensionInfo$MinkowskiPower = MinkowskiPower

  EnsembleClustering = list()
  EnsembleClustering$ClusterRes = Final_Clustering$ClusterRes
  EnsembleClustering$DistMat = Full_Consensus$DistMat
  EnsembleClustering$ConsensusMatrix = Full_Consensus$ConsensusMatrix
  EnsembleClustering$DissimilarityMatrix = Full_Consensus$DissimilarityMatrix
  EnsembleClustering$Embed = Final_MDS_Points
  EnsembleClustering$EmbedDist = Final_Clustering$EmbedDist
  EnsembleClustering$MDS = Final_MDS
  EnsembleClustering$DimensionInfo = DimensionInfo
  EnsembleClustering$HClustRes = Final_Clustering$HClustRes

  return(EnsembleClustering)
}

#Create ensemble clustering from ensemble matrix via MDS + k-means
#' Ensemble clustering by Multidimensional Scaling followed by k-means
#' @description
#' RESTORED FROM EINS v7_1 (EINS_Utils_v6.R). v8 replaced this with
#' `Ensemble_Cluster_MDS_HC()`, which uses `as.dist(1 - consensus)` instead
#' of the Euclidean distance between rows of the dissimilarity matrix, selects
#' the embedding dimension by a stability search rather than setting it equal
#' to the requested cluster count, and clusters the embedding hierarchically
#' rather than with k-means. This original version is kept so the two can be
#' compared on identical data.
#' @param MethodClusterResults A list of dataframes with cluster assignments per
#' multi-omics integration method.
#' @param Sample_IDs Vector of sample names.
#' @param Distance Distance metric to be used for the calculation of the
#' feature distance matrices. Must be one of: `"euclidean"` (default),
#' `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"` or `"minkowski"`.
#' If `"minkowski"` is selected, argument `"MinkowskiPower"` needs to be
#' included.
#' @param MinkowskiPower Power of the Minkowski distance. Default NULL.
#' @param AddConstant Logical, indicating if an additive constant c* should
#' be computed, and added to the non-diagonal dissimilarities such that the
#' modified dissimilarities are Euclidean. Passed to `stats::cmdscale(add = )`.
#' Default TRUE.
#' @param Clusters Number of clusters. Also used as the MDS dimension, as in v7_1.
#' @param KmeansNStart Random starting configurations for k-means. v7_1 used R's
#' default of 1, which made the result non-reproducible; default here is 25.
#' @param Seed Integer seed set immediately before the k-means call, so the
#' result is reproducible in isolation rather than depending on the state of
#' the global random number stream. Default is `NULL`.
#' @returns A list of results. Sample Clustering results, stored in $ClusterRes.
#' Sample distance matrix, stored in $DistMat. MDS final embedding, stored in
#' $Embed.
#' @keywords internal
Ensemble_Cluster_MDS = function(MethodClusterResults,
                                Sample_IDs,
                                Distance = "euclidean",
                                MinkowskiPower = NULL,
                                AddConstant = TRUE,
                                Clusters = NULL,
                                KmeansNStart = 25,
                                Seed = NULL){
  Ens_Sim_Matrix = Ensemble_CoAssociation_Matrix(MethodClusterResults = MethodClusterResults,
                                                 Sample_IDs = Sample_IDs)
  Ens_Dissim_Matrix = 1 - Ens_Sim_Matrix
  if(Distance == "minkowski"){
    if(is.null(MinkowskiPower)){
      Dist_Matrix = stats::dist(Ens_Dissim_Matrix,
                                method = Distance)
    }else{
      Dist_Matrix = stats::dist(Ens_Dissim_Matrix,
                                method = Distance,
                                p = MinkowskiPower)
    }
  }else{
    Dist_Matrix = stats::dist(Ens_Dissim_Matrix,
                              method = Distance)
  }
  MDS = stats::cmdscale(d = Dist_Matrix, k = Clusters, add = AddConstant)
  if(!is.null(Seed)) set.seed(Seed)
  Ens_Clust_Result = stats::kmeans(x = MDS$points, centers = Clusters, nstart = KmeansNStart)[[1]]

  EnsembleClustering = list()
  EnsembleClustering$ClusterRes = Ens_Clust_Result
  EnsembleClustering$DistMat = Dist_Matrix
  EnsembleClustering$Embed = MDS$points

  return(EnsembleClustering)
}

#Median partition of the contributing clusterings
#' Create an ensemble clustering assignment by direct search for the median
#' partition of the contributing clusterings
#' @description
    #' Builds no consensus matrix. Searches instead for the single partition with
    #' the highest total Adjusted Rand Index against every contributing clustering,
    #' which is the median partition of the ensemble under that criterion.
    #' Summarising the contributing clusterings into a co-association matrix
    #' produced by `Ensemble_CoAssociation_Matrix()`, as `Ensemble_Cluster_CHC()`
    #' does, discards information; optimising the agreement criterion directly
    #' does not. The search is a greedy local search over single-sample
    #' reassignments, started from the `Ensemble_Cluster_CHC()` partition and
    #' repeated until no move improves the total or `MaxSweeps` is reached. It is
    #' fully deterministic: the starting partition is deterministic and the sweep
    #' order is fixed, so no seed is needed.
#' @param MethodClusterResults A list of dataframes with cluster assignments per
#' multi-omics integration method, as created by `Ensemble_Cluster_Data()`.
#' @param Sample_IDs Vector of sample names.
#' @param Clusters Number of clusters to create.
#' @param Distance Distance metric used for the starting partition. Must be one
#' of: `"euclidean"` (default), `"maximum"`, `"manhattan"`, `"canberra"`,
#' `"binary"` or `"minkowski"`. If `"minkowski"` is selected, argument
#' `"MinkowskiPower"` needs to be included.
#' @param MinkowskiPower Power of the Minkowski distance. Default is `NULL`.
#' @param Linkage Agglomeration method used for the starting partition. Must be
#' one of: `"ward.D"`, `"ward.D2"` (default), `"single"`, `"complete"`,
#' `"average"`, `"mcquitty"`, `"median"` or `"centroid"`.
#' @param MaxSweeps Maximum number of passes over all samples. Default is 2.
#' @returns A list of results. Sample cluster assignment is stored in $ClusterRes,
#' the co-association matrix is stored in $ConsensusMatrix, the starting sample
#' partition is stored in $StartPartition, the achieved total Adjusted Rand Index
#' is stored in $TotalARI and the number of sweeps performed is stored in $Sweeps.
#' @keywords internal
#' @references Vega-Pons S, Ruiz-Shulcloper J. A survey of clustering ensemble
#' algorithms. International Journal of Pattern Recognition and Artificial
#' Intelligence. 2011;25(3):337-372.
Ensemble_Cluster_MPA = function(MethodClusterResults,
                                Sample_IDs,
                                Clusters = NULL,
                                Distance = "euclidean",
                                MinkowskiPower = NULL,
                                Linkage = "ward.D2",
                                MaxSweeps = 2){
  Ens_Sim_Matrix = Ensemble_CoAssociation_Matrix(MethodClusterResults = MethodClusterResults,
                                                 Sample_IDs = Sample_IDs)
  n_samples = nrow(Ens_Sim_Matrix)

  if(Distance == "minkowski"){
    if(is.null(MinkowskiPower)){
      Dist_Matrix = stats::dist(1 - Ens_Sim_Matrix,
                                method = Distance)
    }else{
      Dist_Matrix = stats::dist(1 - Ens_Sim_Matrix,
                                method = Distance,
                                p = MinkowskiPower)
    }
  }else{
    Dist_Matrix = stats::dist(1 - Ens_Sim_Matrix,
                              method = Distance)
  }
  Start_Partition = stats::cutree(fastcluster::hclust(Dist_Matrix, method = Linkage), k = Clusters)

  Pairs = t(utils::combn(n_samples, 2))
  Base_Pair_Indicators = vapply(MethodClusterResults, function(d){
    Label = d$Cluster
    as.numeric(Label[Pairs[, 1]] == Label[Pairs[, 2]])
  }, numeric(nrow(Pairs)))
  if(is.null(dim(Base_Pair_Indicators))){
    Base_Pair_Indicators = matrix(Base_Pair_Indicators, ncol = 1)
  }
  Base_Pair_Counts = colSums(Base_Pair_Indicators)
  Total_Pairs = n_samples * (n_samples - 1) / 2

  Total_ARI = function(Partition){
    Cand = as.numeric(Partition[Pairs[, 1]] == Partition[Pairs[, 2]])
    Cand_Pairs = sum(Cand)
    Both = as.vector(crossprod(Base_Pair_Indicators, Cand))
    Expected = Cand_Pairs * Base_Pair_Counts / Total_Pairs
    Maximum = (Cand_Pairs + Base_Pair_Counts) / 2
    Index = (Both - Expected) / (Maximum - Expected)
    #0/0 arises when both partitions are trivial in the same way; they then agree
    #perfectly, so the index is 1 rather than NaN.
    Degenerate = (Maximum - Expected) == 0
    if(any(Degenerate)){
      Index[Degenerate] = ifelse(Both[Degenerate] == Cand_Pairs &
                                   Both[Degenerate] == Base_Pair_Counts[Degenerate], 1, 0)
    }
    sum(Index)
  }

  Partition = Start_Partition
  Best_Total = Total_ARI(Partition)
  Sweeps_Done = 0L
  for(sweep in seq_len(MaxSweeps)){
    Sweeps_Done = sweep
    Improved = FALSE
    for(i in seq_len(n_samples)){
      for(g in setdiff(seq_len(Clusters), Partition[i])){
        Candidate = Partition
        Candidate[i] = g
        if(length(unique(Candidate)) != Clusters) next
        Candidate_Total = Total_ARI(Candidate)
        if(Candidate_Total > Best_Total + 1e-12){
          Partition = Candidate
          Best_Total = Candidate_Total
          Improved = TRUE
        }
      }
    }
    if(!Improved) break
  }
  names(Partition) = Sample_IDs

  EnsembleClustering = list()
  EnsembleClustering$ClusterRes = Partition
  EnsembleClustering$ConsensusMatrix = Ens_Sim_Matrix
  EnsembleClustering$StartPartition = Start_Partition
  EnsembleClustering$TotalARI = Best_Total
  EnsembleClustering$Sweeps = Sweeps_Done

  return(EnsembleClustering)
}

#Collect the per-method sample representations used by the embedding ensembles
#' Collect the sample representation of every multi-omics integration method
#' @description
#' Assembles the list of per-method sample representations that
#' `run_Ensemble_Sample_CCA()`, `run_Ensemble_Sample_CAT()`,
#' `run_Ensemble_Sample_CAT_UMAP()` and `run_Ensemble_Sample_CCA_UMAP()` all
#' operate on, so that every one of them consumes an identical starting point
#' rather than rebuilding it. Methods that produce a latent embedding contribute
#' that embedding; COCA and SNF produce no embedding, so their cluster assignment
#' is one-hot encoded instead. Every representation is reordered to `Samples`, so
#' the rows of all of them correspond.
#' @param Samples Vector of sample names, defining the row order of every
#' returned matrix.
#' @param CoordData The `$Multi_Omics$CoordData` list from the main EINS R6 object.
#' Contains sample coordinate data in reduced dimensions.
#' @param Fit The `$Multi_Omics$Fit` list from the main EINS R6 object. Contains
#' all information regarding multi-omics integration method fitting.
#' @param ClusterRes The `$Multi_Omics$ClusterRes` list from the main EINS R6
#' object. Contains the clustering assignments for all multi-omics integration
#' methods.
#' @param Factors Number of latent factors or components the embeddings were
#' created with.
#' @param SNFDistance Which distance used to calculate SNF results to be
#' included in the ensemble clustering. This is needed to reduce the impact
#' of SNF on the ensemble clustering assignment. Must be one of:
#' `"euclidean squared"` (default), `"euclidean"`, `"manhattan"`,
#' `"minkowksi 0.25"`, `"minkowski 0.5"`, `"minkowski 3"` or `"minkowski 4"`.
#' @param GAUDIMethod Which clustering method is used to calculate GAUDI
#' results to be included in the ensemble clustering. This is needed to
#' reduce the impact of GAUDI on the ensemble clustering assignment. Must be
#' one of: `"Kmeans"` (default), `"HDBSCAN"` or `"DBSCAN"`.
#' @param Warn Logical, whether to warn for each method with no representation
#' available. Default is `TRUE`.
#' @returns A named list of matrices, consisting of the embedding per multi-omics
#' integration method. Samples are stored in rows, factors in columns.
#' @keywords internal
Ensemble_Embedding_List = function(Samples,
                                   CoordData,
                                   Fit,
                                   ClusterRes,
                                   Factors,
                                   SNFDistance = "euclidean squared",
                                   GAUDIMethod = "Kmeans",
                                   Warn = TRUE){
  FactorNum = paste0("Factors_", Factors)
  ClustNum = paste0("Clusters_", Factors)
  EmbedList = list()

  #methods with a stored latent embedding
  Coord_Methods = names(CoordData)
  for(method in Coord_Methods){
    if(method == "GAUDI"){
      Embed = CoordData[[method]][[GAUDIMethod]][[FactorNum]]
    }else{
      Embed = CoordData[[method]][[FactorNum]]
    }
    if(!is.null(Embed)){
      Embed = Embed[match(Samples, rownames(Embed)), , drop = FALSE]
      colnames(Embed) = paste0(method, "_Factor_", 1:ncol(Embed))
      EmbedList[[method]] = as.matrix(Embed)
    }else if(Warn){
      warning(paste0("No embedding from ", method, " with ", Factors, " factors available, rerun ", method, " to include this method"))
    }
  }

  #jNMF and iNMF keep their sample loadings in the fit rather than in CoordData
  NMFMethod = c("jNMF", "iNMF")
  for(method in NMFMethod){
    FitW = Fit[[method]][[ClustNum]]$W
    if(!is.null(dim(FitW))){
      FitW = FitW[match(Samples, rownames(FitW)), , drop = FALSE]
      colnames(FitW) = paste0(method, "_Factor_", 1:ncol(FitW))
      EmbedList[[method]] = as.matrix(FitW)
    }else if(Warn){
      warning(paste0("No embedding from ", method, " with ", Factors, " factors available, rerun ", method, " to include this method"))
    }
  }

  #COCA and SNF have no embedding, so the cluster assignment is one-hot encoded
  OneHotMethod = c("COCA", "SNF")
  for(method in OneHotMethod){
    if(method == "SNF"){
      MethodClustRes = ClusterRes[[method]][[SNFDistance]][[ClustNum]]
    }else{
      MethodClustRes = ClusterRes[[method]][[ClustNum]]
    }
    if(!is.null(MethodClustRes)){
      MethodClustRes = MethodClustRes[match(Samples, rownames(MethodClustRes)), , drop = FALSE]
      ClusterLabel = as.numeric(MethodClustRes$Cluster)

      n_samples = length(ClusterLabel)
      #The width must be the number of CLUSTERS, not the number of latent
      #factors. Using Factors made every sample in cluster 10 or above an
      #all-zero row (encoded as belonging to no cluster) whenever
      #Clusters > Factors, and left (Factors - Clusters) empty columns when
      #Clusters < Factors. Fixed in v81.
      n_clusters = max(ClusterLabel, na.rm = TRUE)
      OneHot = matrix(0, nrow = n_samples, ncol = n_clusters)
      for(i in 1:n_samples){
        if(!is.na(ClusterLabel[i])){
          OneHot[i, ClusterLabel[i]] = 1
        }
      }
      rownames(OneHot) = Samples
      colnames(OneHot) = paste0(method, "_Cluster_", 1:n_clusters)
      EmbedList[[method]] = as.matrix(OneHot)
    }else if(Warn){
      warning(paste0("No embedding from ", method, " with ", Factors, " factors available, rerun ", method, " to include this method"))
    }
  }

  return(EmbedList)
}


#Canonical correlation alignment of the per-method representations
#' Align the per-method sample representations by canonical correlation
#' @description
#' Puts every representation into a common space by aligning it to one chosen
#' reference method through canonical correlation analysis, then concatenates the
#' aligned matrices. Without this the representations cannot be concatenated
#' meaningfully, because each method places its own axes arbitrarily.
#' A representation that cannot be aligned is carried through unaligned: the
#' one-hot encodings of COCA and SNF are rank deficient, so `stats::cancor()`
#' returns a coefficient matrix with fewer rows than columns requested and the
#' multiplication fails. That fallback is deliberate and matches the published
#' behaviour of the method.
#' @param EmbedList A named list of per-method embeddings, as created by
#' `Ensemble_Embedding_List()`.
#' @param Reference Name of the method in `EmbedList` to align the others to.
#' @param Samples Vector of sample names, used for the row names of the result.
#' @returns A matrix of the concatenated aligned representations, samples in rows,
#' CCA factor embedding in columns..
#' @keywords internal
Ensemble_Align_CCA = function(EmbedList,
                              Reference,
                              Samples){
  if(Reference %in% names(EmbedList)){
    ref_id = which(names(EmbedList) == Reference)
    EmbedList = c(EmbedList[ref_id], EmbedList[-ref_id])
  }else{
    stop(paste0("No ", Reference, " embeddings available, run ", Reference, " or change Reference variable"))
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

  return(Aligned_Embed)
}


#Standardise and concatenate the per-method representations
#' Standardise and concatenate the per-method sample representations
#' @description
#' Concatenates the per-method representations after standardising each column
#' to zero mean and unit variance, with no alignment step. Standardising is
#' what makes the concatenation meaningful without alignment: methods produce
#' embeddings on very different scales, and an unstandardised concatenation
#' would be dominated by whichever method happens to have the largest variance.
#' This is the control for `Ensemble_Align_CCA()`: the difference between the
#' two isolates what the canonical correlation alignment contributes, since
#' they differ in nothing else. Columns that are constant, and therefore not
#' finite after standardising, are dropped rather than allowed to propagate
#' NaN into the clustering.
#' @param EmbedList A named list of per-method representations, as created by
#' `Ensemble_Embedding_List()`.
#' @param Samples Vector of sample names, used for the row names of the result.
#' @returns A matrix of the concatenated standardised representations, samples in
#' rows.
#' @keywords internal
Ensemble_Concat_Embeddings = function(EmbedList,
                                      Samples){
  Standardise = function(M){
    M = scale(as.matrix(M))
    Finite_Cols = apply(M, 2, function(column) all(is.finite(column)))
    M[, Finite_Cols, drop = FALSE]
  }
  Concat_Embed = do.call(cbind, lapply(EmbedList, Standardise))
  rownames(Concat_Embed) = Samples
  colnames(Concat_Embed) = paste0("CAT_", 1:ncol(Concat_Embed))
  return(Concat_Embed)
}


#Non-linear reduction of a concatenated representation, then k-means
#' Reduce a sample representation with UMAP and cluster the result
#' @description
    #' Applies Uniform Manifold Approximation and Projection to a sample
    #' representation and clusters the reduced coordinates with k-means. Every
    #' other ensemble method in this package combines its inputs linearly - by
    #' averaging co-association counts, by classical multidimensional scaling,
    #' or by canonical correlation - so none of them can recover structure that
    #' is not linearly separable. UMAP can, which is why it is offered as an
    #' alternative reduction step rather than as a replacement for the clustering.
#' @param Embed Matrix of sample coordinates as created by
#' `Ensemble_Concat_Embeddings()`, samples in rows.
#' @param Clusters Number of clusters to create.
#' @param Samples Vector of sample names.
#' @param UMAPNeighbors Size of the local neighborhood used for manifold
#' approximation. Default is 15, reduced automatically if the sample count is too
#' small to support it.
#' @param UMAPDimensions Number of dimensions to reduce to. Default is `NULL`,
#' which uses `Clusters` dimensions, bounded below by 2 and above by the number
#' of available columns.
#' @param UMAPMinDist Effective minimum distance between embedded points. Default
#' is 0.1.
#' @param KmeansNStart Number of random starting configurations for k-means.
#' Default is 25.
#' @param Seed Integer seed set before both the UMAP and the k-means call, so the
#' result is reproducible in isolation. UMAP is stochastic through its stochastic
#' gradient descent, so without this the result depends on the state of the
#' global random number stream. Default is `NULL`, which does not set a seed.
#' @returns A list of results. Cluster assignment in $ClusterRes and the UMAP
#' coordinates in $Embed. Returns `NULL` if UMAP fails, which is possible when the
#' sample count is very small.
#' @keywords internal
#' @references McInnes L, Healy J, Melville J. UMAP: Uniform Manifold
#' Approximation and Projection for Dimension Reduction. arXiv:1802.03426. 2018.
Ensemble_Cluster_UMAP = function(Embed,
                                 Clusters,
                                 Samples,
                                 UMAPNeighbors = 15,
                                 UMAPDimensions = NULL,
                                 UMAPMinDist = 0.1,
                                 KmeansNStart = 25,
                                 Seed = NULL){
  Embed = as.matrix(Embed)
  n_samples = nrow(Embed)
  Neighbors = min(UMAPNeighbors, n_samples - 1)
  if(is.null(UMAPDimensions)){
    UMAPDimensions = max(2, min(Clusters, ncol(Embed)))
  }

  if(!is.null(Seed)) set.seed(Seed)
  UMAP_Embed = try(uwot::umap(Embed,
                              n_neighbors = Neighbors,
                              n_components = UMAPDimensions,
                              min_dist = UMAPMinDist,
                              verbose = FALSE,
                              n_threads = 1,
                              n_sgd_threads = 1),
                   silent = TRUE)
  if(inherits(UMAP_Embed, "try-error")){
    warning(paste0("UMAP failed with ", n_samples, " samples and ", Clusters, " clusters, no result stored"))
    return(NULL)
  }

  UMAP_Embed = as.matrix(UMAP_Embed)
  rownames(UMAP_Embed) = Samples
  colnames(UMAP_Embed) = paste0("UMAP_", 1:ncol(UMAP_Embed))

  if(!is.null(Seed)) set.seed(Seed)
  UMAP_Clust = stats::kmeans(UMAP_Embed, centers = Clusters, nstart = KmeansNStart)$cluster
  names(UMAP_Clust) = Samples

  EnsembleClustering = list()
  EnsembleClustering$ClusterRes = UMAP_Clust
  EnsembleClustering$Embed = UMAP_Embed

  return(EnsembleClustering)
}


#Sample names in the order every ensemble method uses
#' Sample names in the order used by the ensemble methods
#' @description
#' Returns the sample names taken from the most processed data available,
#' falling back through matched, preprocessed, sample and raw data in that
#' order. Every ensemble method resolves its sample order the same way through
#' this function, so their results are row-comparable by construction rather
#' than by coincidence.
#' @param Preprocessed_Omics The `$Preprocessed_Omics` list, which contains the
#' individual omics datasets after preprocessing has been performed.
#' @param Omics The `$Omics` list, which contains the unaltered input omics
#' datasets.
#' @returns A character vector of sample names, ordered as in the first input
#' omics file.
#' @keywords internal
Ensemble_Sample_Names = function(Preprocessed_Omics,
                                 Omics){
  if(!is.null(Preprocessed_Omics$Matched_Data)){
    Data = Preprocessed_Omics$Matched_Data
  } else if(!is.null(Preprocessed_Omics$Data)){
    Data = Preprocessed_Omics$Data
  } else if(!is.null(Omics$Sample_Data)){
    Data = Omics$Sample_Data
  }else{
    Data = Omics$Raw_Data
  }
  return(colnames(Data[[1]]))
}


#Select the contributing base clusterings for the ensemble methods
#' Select the multi-omics cluster assignments that contribute to an ensemble
#' @description
#' Collects the cluster assignments of every multi-omics integration method,
#' keeping only one of the SNF distances and one of the GAUDI clustering
#' algorithms. SNF can be run at several distances, so including all of them
#' would let SNF cast several times as many votes as any other method and
#' dominate the consensus. GAUDI can be run with several clustering
#' algorithms, so including all of them would let GAUDI cast several
#' times as many votes as any other method, and dominate consensus. For
#' @param ClusterRes The `$Multi_Omics$ClusterRes` list. Contains the cluster
#' assignments for all multi-omics integration methods.
#' @param SNFDistance Which distance used to calculate SNF results to be
#' included. Must be one of: `"euclidean squared"` (default), `"euclidean"`,
#' `"manhattan"`, `"minkowksi 0.25"`, `"minkowski 0.5"`, `"minkowski 3"` or
#' `"minkowski 4"`.
#' @param GAUDIMethod Which GAUDI clustering method to take results from.
#' Default is `"Kmeans"`.
#' @returns A list of cluster assignments per method, ready for
#' `Ensemble_Cluster_Data()`.
#' @keywords internal
Ensemble_Base_Clusterings = function(ClusterRes,
                                     SNFDistance = "euclidean squared",
                                     GAUDIMethod = "Kmeans"){
  ClusterResults = list()
  for(method in names(ClusterRes)){
    if(method == "SNF"){
      for(distance in names(ClusterRes$SNF)){
        if(distance == SNFDistance){
          ClusterResults$SNF[[distance]] = ClusterRes$SNF[[distance]]
        }
      }
    }else if(method == "GAUDI"){
      for(clustmeth in names(ClusterRes$GAUDI)){
        if(clustmeth == GAUDIMethod){
          ClusterResults$GAUDI[[clustmeth]] = ClusterRes$GAUDI[[clustmeth]]
        }
      }
    }else{
      ClusterResults[[method]] = ClusterRes[[method]]
    }
  }
  return(ClusterResults)
}

#Hierarchical clustering of scaled factor feature weights using average distance per integration method
#' Ensemble feature clustering using average feature distance between
#' integration methods
#' @description
#' Create an ensemble feature clustering result by combining the feature weights
#' from the multi-omics integration methods. Scaled feature weight matrices of
#' different omics types as calculated by `Feature_Weight_Scaling` are combined
#' per multi-omics integration method. Feature distances are calculated per
#' method, and the average feature distances are calculated between the methods.
#' Hierarchical clustering is performed on the average feature distance matrix.
#' @param ScaledFeaturesData List of scaled feature weight dataframes for the
#' different multi-omics integration methods, by number of factors/components
#' and omics type, as calculated with `Feature_Weight_Scaling()`.
#' @param nFactors Number of factors (components) for which the hierarchical
#' clustering is to be performed. All multi-omics integration methods will need
#' to have been performed with this number of factors/components.
#' @param Distance Distance metric to be used for the calculation of the
#' feature distance matrices. Must be one of: `"euclidean"` (default),
#' `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"` or `"minkowski"`.
#' If `"minkowski"` is selected, argument `"MinkowskiPower"` needs to be
#' included.
#' @param MinkowskiPower Power of the Minkowski distance. Default is `NULL`.
#' @param Linkage Agglomeration method to be used for the hierarchical
#' clustering. Must be one of: `"ward.D"`, `"ward.D2"` (default), `"single"`,
#' `"complete"`, `"average"`, `"mcquitty"`, `"median"` or `"centroid"`.
#' @returns A list of results. Average feature distance matrix, stored in
#' $DistMat. Hierarchical clustering tree, stored in $HClust. Feature dendrogram,
#' stored in $Dendrogram.
#' @keywords internal
Feature_Method_Av_HClust = function(ScaledFeaturesData,
                                    nFactors = NULL,
                                    Distance = "euclidean",
                                    MinkowskiPower = NULL,
                                    Linkage = "ward.D2"
){

  NamedFact = list()
  MethodComb = list()
  OmicsComb = list()
  MethodFeatWeight = list()
  FeatMethods = names(ScaledFeaturesData)
  for(method in FeatMethods){
    FeatFactors = names(ScaledFeaturesData[[method]])
    for(factor in FeatFactors){
      n = as.numeric(strsplit(factor, split = "_")[[1]][2])
      if(n == nFactors){
        FeatOmics = names(ScaledFeaturesData[[method]][[factor]])
        for(omics in FeatOmics){
          Data = subset(ScaledFeaturesData[[method]][[factor]][[omics]], select = 1:(n+1))
          for(i in 1:n){
            Name = paste0(method, "_F", i)
            colnames(Data)[i+1] = Name
          }
          rownames(Data) = NULL
          NamedFact[[method]][[factor]][[omics]] = Data
          MethodComb[[method]][[factor]] = data.table::rbindlist(NamedFact[[method]][[factor]], fill = TRUE)
          OmicsComb[[method]] = purrr::reduce(MethodComb[[method]], dplyr::full_join, by = "Features")
          MethodFeatWeight[[method]] = tibble::column_to_rownames(OmicsComb[[method]], var = "Features")
          MethodFeatWeight[[method]] = t(as.matrix(MethodFeatWeight[[method]]))
        }
      }
    }
  }

  DistList = list()
  for(method in FeatMethods){
    Data = MethodFeatWeight[[method]]
    tData = t(Data)
    #dist creates distance matrix for rows
    if(Distance == "minkowski"){
      if(is.null(MinkowskiPower)){
        DistList[[method]] = stats::dist(tData,
                                         method = Distance)
      }else{
        DistList[[method]] = stats::dist(tData,
                                         method = Distance,
                                         p = MinkowskiPower)
      }
    }else{
      DistList[[method]] = stats::dist(tData,
                                       method = Distance)
    }
  }

  #Average of individual method weigths
  AvFeatDist = Reduce("+", DistList)/length(DistList)

  #hclust
  FeatHClust = fastcluster::hclust(d = AvFeatDist,
                                   method = Linkage)

  FeatDendro = stats::as.dendrogram(FeatHClust)
  FeatDendro |>
    plot()


  HClustRes = list()
  HClustRes$DistMat = AvFeatDist
  HClustRes$HClust = FeatHClust
  HClustRes$Dendrogram = FeatDendro

  return(HClustRes)
}

#Full factor feature dataframe with scaled features
#' Concatenate feature weight results of all multi-omics integration methods
#' into a single matrix.
#' @description
#' Create a single matrix by concatenating the scaled feature weight results as
#' calculated in `Feature_Weight_Scaling()`, with features in rows and factors
#' (components) per multi-omics integration method in columns. To be used in
#' `Feature_Concat_HClust()`.
#' @param ScaledFeaturesData List of scaled feature weight dataframes for the
#' different multi-omics integration methods, by number of factors/components
#' and omics type, as calculated with `Feature_Weight_Scaling()`.
#' @param nFactors Number of factors (components) for which the hierarchical
#' clustering is to be performed. All multi-omics integration methods will need
#' to have been performed with this number of factors/components.
#' @returns A single matrix of scaled feature weights. Features in rows, factors
#' or components per multi-omics integration method in columns.
#' @keywords internal
Concat_Feature_Matrix = function(ScaledFeaturesData,
                                 nFactors = NULL){

  NamedFact = list()
  MethodComb = list()
  OmicsComb = list()

  FeatMethods = names(ScaledFeaturesData)

  for(method in FeatMethods){
    FeatFactors = names(ScaledFeaturesData[[method]])

    for(factor in FeatFactors){
      n = as.numeric(strsplit(factor, split = "_")[[1]][2])

      if(n == nFactors){
        FeatOmics = names(ScaledFeaturesData[[method]][[factor]])

        for(omics in FeatOmics){
          Data = subset(ScaledFeaturesData[[method]][[factor]][[omics]], select = 1:(n+1))

          for(i in 1:n){
            Name = paste0(method, "_F", i)
            colnames(Data)[i+1] = Name
          }

          rownames(Data) = NULL
          NamedFact[[factor]][[omics]][[method]] = Data
          MethodComb[[factor]][[omics]] = purrr::reduce(NamedFact[[factor]][[omics]], dplyr::full_join, by = "Features")
          OmicsComb[[factor]] = data.table::rbindlist(MethodComb[[factor]], fill = TRUE)
          FullFeatWeight = purrr::reduce(OmicsComb, dplyr::full_join, by = "Features")
          FullFeatWeight = tibble::column_to_rownames(FullFeatWeight, var = "Features")
        }
      }
    }
  }

  FullFeatWeight = t(as.matrix(FullFeatWeight))

  return(FullFeatWeight)
}

#Hierarchical clustering of all factor concatenated feature weights
#' Ensemble feature and factor clustering using concatenated feature weight
#' matrix
#' @description
#' Create an ensemble feature clustering result by combining the feature weights
#' from the multi-omics integration methods. Single matrix of scaled feature
#' weights from all multi-omics integration methods as calculated by
#' `Concat_Feature_Matrix()` is needed. Feature and factor distances are
#' calculated on this matrix. Hierarchical clustering is performed for the
#' feature and factor distances.
#' @param FullFeatureData Concatenated matrix of scaled feature weights from all
#' multi-omics integration methods as calculated by `Concat_Feature_Matrix()`.
#' @param Distance Distance metric to be used for the calculation of the
#' feature distance matrices. Must be one of: `"euclidean"` (default),
#' `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"` or `"minkowski"`.
#' If `"minkowski"` is selected, argument `"MinkowskiPower"` needs to be
#' included.
#' @param MinkowskiPower Power of the Minkowski distance. Default is `NULL`.
#' @param Linkage Agglomeration method to be used for the hierarchical
#' clustering. Must be one of: `"ward.D"`, `"ward.D2"` (default), `"single"`,
#' `"complete"`, `"average"`, `"mcquitty"`, `"median"` or `"centroid"`.
#' @returns A list of results. Concatenated feature distance matrix, stored in
#' $FeatureDist. Concatenated factor distance matrix, stored in $FactorDist.
#' Feature hierarchical clustering tree, stored in $FeatureHClust. Factor
#' hierarchical clustering tree, stored in $FactorHClust. Feature dendrogram,
#' stored in $FeatureDendro. Factor dendrogram, stored in $FactorDendro.
#' @keywords internal
Feature_Concat_HClust = function(FullFeatureData,
                                 Distance = "euclidean",
                                 MinkowskiPower = NULL,
                                 Linkage = "ward.D2"
){
  HClustRes = list()
  Data = FullFeatureData
  tData = t(Data)
  #dist creates distance matrix for rows
  if(Distance == "minkowski"){
    if(is.null(MinkowskiPower)){
      FeatDist = stats::dist(tData,
                             method = Distance)
      FactDist = stats::dist(Data,
                             method = Distance)
    }else{
      FeatDist = stats::dist(tData,
                             method = Distance,
                             p = MinkowskiPower)
      FactDist = stats::dist(Data,
                             method = Distance,
                             p = MinkowskiPower)
    }
  }else{
    FeatDist = stats::dist(tData,
                           method = Distance)
    FactDist = stats::dist(Data,
                           method = Distance)
  }

  #hclust
  FeatHClust = fastcluster::hclust(d = FeatDist,
                                   method = Linkage)

  FactHClust = fastcluster::hclust(d = FactDist,
                                   method = Linkage)

  #plot dendrograms
  FeatDendro = stats::as.dendrogram(FeatHClust)
  FeatDendro |>
    plot()
  FactDendro = stats::as.dendrogram(FactHClust)
  FactDendro |>
    plot()

  HClustRes$FactorDist = FactDist
  HClustRes$FeatureDist = FeatDist
  HClustRes$FeatureHClust = FeatHClust
  HClustRes$FactorHClust = FactHClust
  HClustRes$FeatureDendro = FeatDendro
  HClustRes$FactorDendro = FactDendro
  return(HClustRes)
}

#Optimal cluster number for feature distance matrix
#' Calculate the optimal cluster number for features or factors using the
#' hierarchical clustering result.
#' @description
#' Calculate the optimal number of features or factors using the hierarchical
#' clustering trees and distance matrices calculated in either
#' `Feature_Method_Av_HClust()` or `Feature_Concat_HClust()`. Uses the
#' as.clustrange function from the WeightedCluster package.
#' @param HClust Feature or factor hierarchical clustering tree as calculated in
#' either `Feature_Method_Av_HClust()` or `Feature_Concat_HClust()`.
#' @param Dist Feature or factor distance matrix as calculated in
#' `Feature_Method_Av_HClust()` or `Feature_Concat_HClust()`.
#' @param WeightedClusterStat The statistics included to calculate the
#' optimal cluster number. Must be one of: `"all"` (default), `"noCH"` (all
#' statistics except `"CH"` and `"CHsq"`), `"PBC"`, `"HG"`, `"HGSD"`,
#' `"ASW"`, `"ASWw"`, `"CH"`, `"R2"`, `"CHsq"`, `"R2sq"` or `"HC"`.
#' @param MaxClusters The maximum number of feature or factor clusters to be
#' considered by the statistics.
#' @returns The optimal number of feature or factor clusters.
#' @keywords internal
#' @references Studer M (2013). “WeightedCluster Library Manual: A practical
#' guide to creating typologies of trajectories in the social sciences with R.”
#' LIVES Working Papers 24.
Feature_Optimal_Clusters = function(HClust,
                                    Dist,
                                    WeightedClusterStat,
                                    MaxClusters){
  WClust = WeightedCluster::as.clustrange(object = HClust,
                                          diss = Dist,
                                          weights = NULL,
                                          R = 1,
                                          ncluster = MaxClusters,
                                          stat = WeightedClusterStat)
  WClustSum = summary(WClust)
  K = round(median(WClustSum$`1. N groups`), digits = 0)
  OptClustFeat = list()
  OptClustFeat$OptClust$Fit = WClust
  OptClustFeat$OptClust$K = K
  return(OptClustFeat)
}

#feature hierarchical clustering dendrogram grouping
#' Cut feature or factor dendrograms into clusters
#' @description
#' Cut feature of factor dendrograms as created in either
#' `Feature_Method_Av_HClust()` or `Feature_Concat_HClust()` to create cluster
#' assignments.
#' @param FeatHClust Hierarchical clustering results for features or factors
#' as calculated in either `Feature_Method_Av_HClust()` or
#' `Feature_Concat_HClust()`.
#' @param Clusters Number of clusters in which to cut the feature or factor
#' hierarchical clustering tree.
#' @returns A list of results. Feature or factor cluster assignment, stored in
#' $Clusters. Feature or factor dendrogram with branches colored for the number
#' of clusters, stored in $Dendro.
#' @keywords internal
Feature_Dendro_Clustering = function(FeatHClust,
                                     Clusters){
  #cuttree
  FeatClust = dendextend::cutree(tree = FeatHClust,
                                 k = Clusters,
                                 order_clusters_as_data = FALSE)
  FeatClustRes = data.frame(Cluster = FeatClust,
                            row.names = FeatHClust$labels[FeatHClust$order],
                            stringsAsFactors = FALSE)
  FeatDendro = stats::as.dendrogram(FeatHClust)
  FeatDendroPlot = FeatDendro |>
    dendextend::set("labels_colors", k = Clusters) |>
    dendextend::set("branches_k_color", k = Clusters)

  FeatClust = list()
  FeatClust$Dendro = FeatDendroPlot
  FeatClust$Clusters = FeatClustRes
  return(FeatClust)
}

#Enrichment analysis
#' Feature overrepresentation analysis on the clustered ensemble feature results
#' @description
#' Perform overrepresentation analysis from package ClusterProfiler on the
#' clustered ensembled features to determine if features with similar biological
#' functions group together in across multi-omics integration methods. Uses
#' enrichment GO categories from the genome wide annotation database for human
#' as provided in package org.Hs.eg.db. Performed per cluster and per omics type.
#' Uses feature cluster assignments created in `Feature_Dendro_Clustering()`.
#' @param FullList Vector of all features of the analyzed omics type.
#' @param SelectList Vector of features of the analyzed omics type in the
#' cluster.
#' @param GeneNameType Keytype of the feature name for the omics dataset. Must
#' be one of: `"ENTREZID"`, `"PFAM"`, `"IPI"`, `"PROSITE"`, `"ACCNUM"`,
#' `"ALIAS"`, `"CHR"`, `"CHRLOC"`, `"CHRLOCEND"`, `"ENZYME"`, `"MAP"`, `"PATH"`,
#' `"PMID"`, `"REFSEQ"`, `"SYMBOL"`, `"UNIGENE"`, `"ENSEMBL"`, `"ENSEMBLPROT"`,
#' `"ENSEMBLTRANS"`, `"GENENAME"`, `"UNIPROT"`, `"GO"`, `"EVIDENCE"`,
#' `"ONTOLOGY"`, `"GOALL"`, `"EVIDENCEALL"`, `"ONTOLOGYALL"`, `"OMIM"` or
#' `"UCSCKG"`.
#' @param Subontologies Which ontology to use. Must be one of: `"ALL"` (default),
#' `"MF"`, `"BP"` or `"CC"`.
#' @param pValue Cutoff value for the p-value. Default is 0.05.
#' @param pAdjustment Method for the p-value adjustment. Must be one of: `"BH"`
#' (default), `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BY"`, `"fdr"`
#' or `"none"`.
#' @param qValue Cutoff for the q-value. Default is 0.2.
#' @param MinGenesPerTerm Minimal number of genes annotated per ontology term to
#' be included for testing. Default is 10.
#' @param MaxGenesPerTerm Maximum number of genes annotated per ontology term to
#' be included for testing. Default is 500.
#' @returns A list of results stored in class `"enrichResult"`. Enrichment
#' analysis result, stored in $result. p-value cutoff value, stored in
#' $pvalueCutoff. p-value adjustment method, stored in $pAdjustMethod. q-value
#' cutoff value, stored in $qvalueCutoff. Organism for which enrichment is
#' performed, can only be `"human"`, stored in $organism. Ontology used, stored
#' in $ontology. Feature IDs in cluster, stored in $gene. Feature keytype,
#' stored in $keytype. All feature ID's, stored in $universe. Feature and
#' category association, stored in $geneInCategory. Mapping of feature to
#' Symbol, stored in $gen2Symbol. Feature sets, stored in $geneSets. Logical
#' flag of feature ID in symbol or not, stored in $readable.
#' @keywords internal
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
ORA_GO = function(FullList,
                  SelectList,
                  GeneNameType,
                  Subontologies = "ALL",
                  pValue = 0.05,
                  pAdjustment = "BH",
                  qValue = 0.2,
                  MinGenesPerTerm = 10,
                  MaxGenesPerTerm = 500
){
  ENSFullList = clusterProfiler::bitr(geneID = FullList,
                                      fromType = GeneNameType,
                                      toType = "ENSEMBL",
                                      OrgDb = org.Hs.eg.db)
  ENSSelectList = clusterProfiler::bitr(geneID = SelectList,
                                        fromType = GeneNameType,
                                        toType = "ENSEMBL",
                                        OrgDb = org.Hs.eg.db)

  GOEnriched = clusterProfiler::enrichGO(gene = ENSSelectList$ENSEMBL,
                                         universe = ENSFullList$ENSEMBL,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = "ENSEMBL",
                                         ont = Subontologies,
                                         pvalueCutoff = pValue,
                                         pAdjustMethod = pAdjustment,
                                         qvalueCutoff = qValue,
                                         minGSSize = MinGenesPerTerm,
                                         maxGSSize = MaxGenesPerTerm)
  return(GOEnriched)
}

#----------------Visualization----------------

#Heatmap function with cluster assignment
#' Interactive heatmap with sample clustering result
#' @description
#' Interactive heatmap, with sample clustering result.
#' @param Data Omics matrix with features in rows and samples in columns.
#' @param SampleClusters Dataframe with sample clustering result with samples in
#' rows and cluster assignment in the column.
#' @param FeatureClusters Whether to order the features on the Y-axis based on
#' hierarchical clustering. Default is `TRUE`. If `FALSE`, the features will be
#' ordered as in `Data`.
#' @param YFontSize Font size of the feature labels on the Y-axis. Default is 5.
#' @param XFontSize Font size of the sample labels on the X-axis. Default is 5.
#' @param Title Title for the heatmap.
#' @returns Interactive heatmap with cluster assignment on X-axis.
#' @keywords internal
Heatmap_Cluster = function(Data,
                           SampleClusters = NULL,
                           FeatureClusters = TRUE,
                           YFontSize = 5,
                           XFontSize = 5,
                           Title = NULL) {
  #Set clusters
  SampleClustersVector = as.vector(SampleClusters$Cluster)
  ClusterOrder = order(SampleClustersVector)
  SampleClusters$Cluster = as.factor(SampleClusters$Cluster)


  #reorder Data columns based on clustering
  Data = Data[, ClusterOrder]
  SampleClusters = SampleClusters[ClusterOrder, ]

  #plot data
  Plot = heatmaply::heatmaply(x = Data,
                              plot_method = "ggplot",
                              Colv = FALSE,
                              Rowv = FeatureClusters,
                              col_side_colors = SampleClusters,
                              main = Title,
                              xlab = "Samples",
                              ylab = "Features",
                              fontsize_row = XFontSize,
                              fontsize_col = YFontSize,
                              label_names = c("Feature", "Sample", "Value"))

  Plot = plotly::layout(p = Plot,
                        plot_bgcolor = "rgba(0,0,0,0)",
                        paper_bgcolor = "rgba(0,0,0,0)")
  return(Plot)
}

#Create a single file with sample metadata for the matched samples only
#' Combine sample metadata from different omics for samples matches across omics
#' @description
#' Create a dataframe with identical metadata between different omics metadata
#' files, for the samples which exist across all omics.
#' @param MetadataList List of metadata files sorted by omics.
#' @returns Dataframes with matching metadata in matching samples.
#' @keywords internal
Match_Metadata = function(MetadataList){
  #wanted output: parts of different metadata files which are identical and present for matched samples kept in samples in rows, features in columns dataframe
  tMetadataList = list()
  Omics_Names = names(MetadataList)

  #create list of dataframes with samples in columns
  for(i in Omics_Names){
    tMetadataList[[i]] = data.table::transpose(MetadataList[[i]])
    colnames(tMetadataList[[i]]) = rownames(MetadataList[[i]])
    rownames(tMetadataList[[i]]) = colnames(MetadataList[[i]])
    tMetadataList[[i]] = tibble::rownames_to_column(tMetadataList[[i]], "Features")
  }
  #rowbind all dataframes from the list together
  FullMetaList = dplyr::bind_rows(tMetadataList)

  #Keep only matched samples, no NA after rowbind
  MetaSampleList = list()
  MetaSampleList = lapply(tMetadataList, colnames)
  MetaMatchedSampleList = Reduce(intersect, MetaSampleList)
  MatchMetaList = FullMetaList[MetaMatchedSampleList]
  #keep 1 set of the duplicated columns
  DupMetaList = dplyr::distinct(MatchMetaList[duplicated(MatchMetaList),])
  rownames(DupMetaList) = NULL
  #have Features column be used as rownames
  DupMetaList = tibble::column_to_rownames(DupMetaList, "Features")
  #retranspose the data with samples in rows
  Matched_Metadata_List = data.table::transpose(DupMetaList)
  colnames(Matched_Metadata_List) = rownames(DupMetaList)
  rownames(Matched_Metadata_List) = colnames(DupMetaList)
  return(Matched_Metadata_List)
}

#multi omics heatmap without cluster annotations
#' Non-interactive heatmap with all omics in single plot
#' @description
#' A non-interactive heatmap with all omics in a single plot, and sample
#' metadata features plotted on the X-axis.
#' @param OmicsData List of omics matrices, with features in rows and samples
#' in columns.
#' @param LegendNames  Names for the legends of the different omics heatmaps.
#' If `NULL`, the omics names as used in the EINS object are used. Otherwise,
#' a string vector with the same length as `OmicsData` should be provided.
#' @param MetadataColumn Names of the sample metadata features to be displayed
#' on the heatmaps. These metadata features need to be available in the metadata
#' files for all omics types. A string vector with the names of the metadata
#' features as in `MetaData` should be provided.
#' @param MetaData List of metadata dataframes for all omics.
#' @returns Multi-omics heatmap with metadata features on the X-axis.
#' @keywords internal
MultiOmicsHeatmap = function(OmicsData,
                             LegendNames = NULL, #vector with legend names per omics (ordered as omics in OmicsData)
                             MetadataColumn = NULL,
                             MetaData = NULL
){
  if(!is.null(MetadataColumn)){
    MatchedMetaData = Match_Metadata(MetaData)
    ReqMetaData = MatchedMetaData[MetadataColumn]
    MetaNames = colnames(ReqMetaData)
    ReqMetaData = tibble::rownames_to_column(ReqMetaData, var = "Sample")

    Meta_Short = list(Meta1 = RColorBrewer::brewer.pal(9, "Set1"),
                      Meta2 = RColorBrewer::brewer.pal(8, "Set2"),
                      Meta3 = RColorBrewer::brewer.pal(8, "Dark2"),
                      Meta4 = RColorBrewer::brewer.pal(8, "Accent"),
                      Meta5 = RColorBrewer::brewer.pal(12, "Set3"),
                      Meta6 = RColorBrewer::brewer.pal(9, "Pastel1"))
    Meta_Long = list(Meta1 = Polychrome::glasbey.colors(n = 32),
                     Meta2 = Polychrome::palette36.colors(n = 36),
                     Meta3 = Polychrome::dark.colors(n = 24),
                     Meta4 = Polychrome::alphabet.colors(n = 26),
                     Meta5 = Polychrome::sky.colors(n = 24),
                     Meta6 = Polychrome::green.armytage.colors(n = 26))

    MetadataColorList =  list()

    for(i in 1:length(MetadataColumn)){
      Feature = MetadataColumn[i]
      Meta_Factor = as.factor(ReqMetaData[, Feature])
      Meta_Name = paste0("Meta", i)
      if(nlevels(Meta_Factor) <= 10){
        MetadataColorList[[Meta_Name]] = Meta_Short[[i]]
      }else{
        MetadataColorList[[Meta_Name]] = Meta_Long[[i]]
      }
    }

    #list of vectors for metadata annotation
    AnnoColor = list()
    for(i in 2:ncol(ReqMetaData)){
      AnnoColor[[i - 1]] = MetadataColorList[[i - 1]][1:nrow(unique(ReqMetaData[i]))]
      names(AnnoColor[[i - 1]]) = as.vector(unique(ReqMetaData[i]))[[1]]
    }
    names(AnnoColor) = MetaNames

    AnnoInfo = ReqMetaData
    AnnoInfo = subset(AnnoInfo, select = -Sample)

    HeatmapAnno = ComplexHeatmap::HeatmapAnnotation(df = AnnoInfo,
                                                    col = AnnoColor,
                                                    border = FALSE)

  } else{
    HeatmapAnno = NULL
  }



  Omics1 = c("#FF0000", "#FFFF00", "#00FF00")
  Omics2 = c("#0000FF", "#808080", "#FFFF00")
  Omics3 = c("#FF8000", "#FFFFFF", "#7502f7")
  Omics4 = c("#FF00FF", "#808080", "#00FF00")
  Omics5 = c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
  Omics6 = c("#00FF00", "#FFFFFF", "#7502f7")
  ColorPalettes = list(Omics1 = Omics1, Omics2 = Omics2, Omics3 = Omics3, Omics4 = Omics4, Omics5 = Omics5, Omics6 = Omics6)

  OmicsNames = names(OmicsData)
  if(is.null(LegendNames)){
    LegendNames = OmicsNames
  }
  OmicsHeatmaps = list()
  for(i in 1:length(OmicsData)){
    RowTitle = paste0(OmicsNames[i])
    OmicsHeatmaps[[i]] = ComplexHeatmap::Heatmap(matrix = OmicsData[[i]],
                                                 row_title = RowTitle,
                                                 name = LegendNames[i],
                                                 cluster_columns = FALSE,
                                                 show_column_dend = FALSE,
                                                 show_column_names = TRUE,
                                                 cluster_rows = TRUE,
                                                 show_row_dend = TRUE,
                                                 show_row_names = FALSE,
                                                 col = grDevices::colorRampPalette(ColorPalettes[[i]])(5),
                                                 top_annotation = switch((i == 1) + 1, NULL, HeatmapAnno),
                                                 height = 5,
                                                 width = 10
    )
  }
  MultiOmicsHeatmap = Reduce(ComplexHeatmap::"%v%", OmicsHeatmaps)
  MultiOmicsItems = list(Heatmap = MultiOmicsHeatmap, Annotation = HeatmapAnno)
  return(MultiOmicsItems)
}

#multi omics heatmap with cluster annotation
#' Non-interactive heatmap with cluster annotation for all omics in single plot
#' @description
#' A non-interactive heatmap with all omics in a single plot, and multi-omics
#' integration or ensemble sample cluster assignment as well as sample metadata
#' features plotted on the X-axis.
#' @param OmicsData List of omics matrices, with features in rows and samples
#' in columns.
#' @param ClusterRes Dataframe with a sample cluster assignment.
#' @param LegendNames Names for the legends of the different omics heatmaps.
#' If `NULL`, the omics names as used in the EINS object are used. Otherwise,
#' a string vector with the same length as `OmicsData` should be provided.
#' @param MetadataColumn Names of the sample metadata features to be displayed
#' on the heatmaps. These metadata features need to be available in the metadata
#' files for all omics types. A string vector with the names of the metadata
#' features as in `MetaData` should be provided.
#' @param MetaData List of metadata dataframes for all omics.
#' @returns Multi-omics heatmap with cluster assignment and metadata features
#' on the X-axis.
#' @keywords internal
MultiOmicsHeatmapClustered = function(OmicsData,
                                      ClusterRes,
                                      LegendNames = NULL,
                                      MetadataColumn = NULL,
                                      MetaData = NULL){
  #create color palettes
  ClusterColors = RColorBrewer::brewer.pal(12, "Paired")

  #vector for cluster annotation
  ClusterRes = tibble::rownames_to_column(ClusterRes, var = "Sample")
  ClusColVec = ClusterColors[1:length(unique(ClusterRes$Cluster))]
  names(ClusColVec) = paste0(unique(ClusterRes$Cluster))

  if(!is.null(MetadataColumn)){
    #create df with only wanted metadata columns
    MatchedMetaData = Match_Metadata(MetaData)
    ReqMetaData = MatchedMetaData[MetadataColumn]
    MetaNames = colnames(ReqMetaData)
    ReqMetaData = tibble::rownames_to_column(ReqMetaData, var = "Sample")

    Meta_Short = list(Meta1 = RColorBrewer::brewer.pal(9, "Set1"),
                      Meta2 = RColorBrewer::brewer.pal(8, "Set2"),
                      Meta3 = RColorBrewer::brewer.pal(8, "Dark2"),
                      Meta4 = RColorBrewer::brewer.pal(8, "Accent"),
                      Meta5 = RColorBrewer::brewer.pal(12, "Set3"),
                      Meta6 = RColorBrewer::brewer.pal(9, "Pastel1"))
    Meta_Long = list(Meta1 = Polychrome::glasbey.colors(n = 32),
                     Meta2 = Polychrome::palette36.colors(n = 36),
                     Meta3 = Polychrome::dark.colors(n = 24),
                     Meta4 = Polychrome::alphabet.colors(n = 26),
                     Meta5 = Polychrome::sky.colors(n = 24),
                     Meta6 = Polychrome::green.armytage.colors(n = 26))

    MetadataColorList =  list()

    for(i in 1:length(MetadataColumn)){
      Feature = MetadataColumn[i]
      Meta_Factor = as.factor(ReqMetaData[, Feature])
      Meta_Name = paste0("Meta", i)
      if(nlevels(Meta_Factor) <= 10){
        MetadataColorList[[Meta_Name]] = Meta_Short[[i]]
      }else{
        MetadataColorList[[Meta_Name]] = Meta_Long[[i]]
      }
    }

    #list of vectors for metadata annotation
    AnnoColor = list()
    for(i in 2:ncol(ReqMetaData)){
      AnnoColor[[i - 1]] = MetadataColorList[[i - 1]][1:nrow(unique(ReqMetaData[i]))]
      names(AnnoColor[[i - 1]]) = as.vector(unique(ReqMetaData[i]))[[1]]
    }
    names(AnnoColor) = MetaNames

    AnnoInfo = merge(ClusterRes, ReqMetaData, by = "Sample", sort = FALSE)
    ClusterOrder = order(ClusterRes$Cluster)
    AnnoInfo = AnnoInfo[ClusterOrder, ]
    AnnoInfo = subset(AnnoInfo, select = -Sample)
    AnnoColor[["Cluster"]] = ClusColVec

  } else{
    AnnoInfo = ClusterRes
    ClusterOrder = order(ClusterRes$Cluster)
    AnnoInfo = AnnoInfo[ClusterOrder, ]
    AnnoInfo = subset(AnnoInfo, select = -Sample)
    AnnoColor = list("Cluster" = ClusColVec)
  }

  HeatmapAnno = ComplexHeatmap::HeatmapAnnotation(df = AnnoInfo,
                                                  col = AnnoColor,
                                                  border = FALSE)

  Omics1 = c("#FF0000", "#FFFF00", "#00FF00")
  Omics2 = c("#0000FF", "#808080", "#FFFF00")
  Omics3 = c("#FF8000", "#FFFFFF", "#7502f7")
  Omics4 = c("#FF00FF", "#808080", "#00FF00")
  Omics5 = c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
  Omics6 = c("#00FF00", "#FFFFFF", "#7502f7")
  ColorPalettes = list(Omics1 = Omics1, Omics2 = Omics2, Omics3 = Omics3, Omics4 = Omics4, Omics5 = Omics5, Omics6 = Omics6)

  OmicsNames = names(OmicsData)
  if(is.null(LegendNames)){
    LegendNames = OmicsNames
  }
  OmicsHeatmaps = list()
  for(i in 1:length(OmicsData)){
    RowTitle = paste0(OmicsNames[i])
    OrderedOmicsData = OmicsData[[i]][,ClusterOrder]
    OmicsHeatmaps[[i]] = ComplexHeatmap::Heatmap(matrix = OrderedOmicsData,
                                                 row_title = RowTitle,
                                                 name = LegendNames[i],
                                                 cluster_columns = FALSE,
                                                 show_column_dend = FALSE,
                                                 show_column_names = TRUE,
                                                 cluster_rows = TRUE,
                                                 show_row_dend = TRUE,
                                                 show_row_names = FALSE,
                                                 col = grDevices::colorRampPalette(ColorPalettes[[i]])(5),
                                                 top_annotation = switch((i == 1) + 1, NULL, HeatmapAnno),
                                                 height = 5,
                                                 width = 10
    )
  }
  MultiOmicsHeatmapAnno = Reduce(ComplexHeatmap::"%v%", OmicsHeatmaps)
  return(MultiOmicsHeatmapAnno)
}

#interactive scatterplot
#' Interactive sample scatterplot
#' @description
#' Interactive scatterplot of samples in the factors (components/dimensions),
#' colored by cluster.
#' @param Data Dataframe with samples in rows, columns with factors
#' (components/dimensions), a column with cluster assignment and an optional
#' column with a sample metadata feature for the shape of the dots.
#' @param XAxis Integer, factor (component/dimension) to plot on the X-axis.
#' Default is 1.
#' @param YAxis Integer, factor (component/dimension) to plot on the Y-axis.
#' Default is 2.
#' @param Width Width of the scatterplot. Default is 5.
#' @param Height Height of the scatterplot. Default is 5.
#' @param YFontSize Size of the Y-axis title. Default is 15.
#' @param XFontSize Size of the X-axis title. Default is 15.
#' @param Title Title of the scatterplot.
#' @returns An interactive sample scatterplot, colored by cluster and,
#' optionally, with shapes determined by a metadata feature.
#' @keywords internal
Scatterplot = function(Data,
                       XAxis = XAxis,
                       YAxis = YAxis,
                       Width = 5,
                       Height = 5,
                       YFontSize = 15,
                       XFontSize = 15,
                       Title){
  Plot = plotly::plot_ly(data = Data,
                         x = Data[,XAxis],
                         y = Data[,YAxis],
                         color = Data$Metadata,
                         colors = "Set1",
                         symbol = Data$Cluster,
                         width = Width,
                         height = Height)

  Plot = plotly::add_trace(p = Plot,
                           type = "scatter",
                           mode = "markers",
                           text = ~paste("Sample: ", rownames(Data))
  )
  Plot = plotly::layout(p = Plot,
                        title = list(text = Title),
                        xaxis = list(title = list(text = names(Data[XAxis]),
                                                  font = list(size = XFontSize))),
                        yaxis = list(title = list(text = names(Data[YAxis]),
                                                  font = list(size = YFontSize))),
                        legend = list(title = list(text = "Cluster"))
  )
  return(Plot)
}

#Create colored table of cluster assignment and metadata
#' Table with cells colored for sample cluster assignment and metadata features
#' @description
#' Table with columns of sample cluster assignment and metadata features, with
#' each column having a different color palette.
#' @param MetadataColumns Names of the sample metadata features to be included
#' in the table. These metadata features need to be available in the metadata
#' files for all omics types. A string vector with the names of the metadata
#' features as in the metadata files should be provided.
#' @param MatchedMetadata A dataframe of matched metadata features across omics,
#' as provided by `Match_Metadata()`.
#' @param ClusterResults Dataframe with sample cluster assignment.
#' @returns A cluster assignment and metadata feature table.
#' @keywords internal
Plot_Clusters_Metadata = function(MetadataColumns = NULL,
                                  MatchedMetadata,
                                  ClusterResults){
  #keep only useful metadata columns
  ReqMetadata = MatchedMetadata[MetadataColumns]
  #merge metadata and clustering results
  ClustMeta = merge(ClusterResults, ReqMetadata, by = "row.names")
  #rename Row.names to sample id
  names(ClustMeta)[names(ClustMeta) == "Row.names"] <- "Sample ID"
  #Order data based on cluster number
  ClustMeta = dplyr::arrange(ClustMeta, Cluster)

  #change class of character columns to factors
  ClustMeta[sapply(ClustMeta, is.character)] = lapply(ClustMeta[sapply(ClustMeta, is.character)], as.factor)
  Meta_Short = list(Greys = RColorBrewer::brewer.pal(9, "Greys"),
                    Clusters = RColorBrewer::brewer.pal(12, "Paired"),
                    Meta1 = RColorBrewer::brewer.pal(9, "Set1"),
                    Meta2 = RColorBrewer::brewer.pal(8, "Set2"),
                    Meta3 = RColorBrewer::brewer.pal(8, "Dark2"),
                    Meta4 = RColorBrewer::brewer.pal(8, "Accent"),
                    Meta5 = RColorBrewer::brewer.pal(12, "Set3"),
                    Meta6 = RColorBrewer::brewer.pal(9, "Pastel1"))
  Meta_Long = list(Greys = RColorBrewer::brewer.pal(9, "Greys"),
                   Clusters = RColorBrewer::brewer.pal(12, "Paired"),
                   Meta1 = Polychrome::glasbey.colors(n = 32),
                   Meta2 = Polychrome::palette36.colors(n = 36),
                   Meta3 = Polychrome::dark.colors(n = 24),
                   Meta4 = Polychrome::alphabet.colors(n = 26),
                   Meta5 = Polychrome::sky.colors(n = 24),
                   Meta6 = Polychrome::green.armytage.colors(n = 26))
  PlotClustMeta = gt::gt(ClustMeta)
  for(i in 2:ncol(ClustMeta)){
    if(nlevels(ClustMeta[[i]]) <= 10){
      PlotClustMeta = gt::data_color(data = PlotClustMeta, columns = i, palette = Meta_Short[[i]])
    }else{
      PlotClustMeta = gt::data_color(data = PlotClustMeta, columns = i, palette = Meta_Long[[i]])
    }
  }
  return(PlotClustMeta)
}

#Data manipulation for the Sankey plot comparing differing cluster number results
#' Create dataframe to be used for Sankey plot comparing cluster numbers
#' @description
#' Combine the cluster results for a single multi-omics method performed for
#' different cluster numbers in a single dataframe. For each sample with each
#' increase in cluster number, the sample is given a source and target value,
#' which can be used to display the flow in the Sankey plot. Each sample also
#' gets a color assigned based on a user-defined metadata feature.
#' @param ClusterResList List of cluster assignment dataframes.
#' @param MetadataColumn String of the sample metadata feature for which the
#' samples should be colored. The metadata feature needs to be available
#' in the metadata files for all omics types.
#' @returns List with results. $Links, dataframe with samples in rows, and
#' columns indicating source cluster and target cluster, as well as metadata
#' color. $Nodes, dataframe with source and target clusters.
#' @keywords internal
Data_Manipulation_Sankey_Clusters = function(ClusterResList,
                                             MetadataColumn){
  #identify all cluster numbers used in ClusterResList
  Sankey_Clustnum = names(ClusterResList)
  #Rename the Cluster column in the ClusterResList to reflect the cluster number
  for(clustnum in Sankey_Clustnum){
    names(ClusterResList[[clustnum]])[names(ClusterResList[[clustnum]]) == "Cluster"] <- clustnum
  }
  #combine all ClusterRes dataframes by column
  Sankey_Full_Clust_df = do.call(cbind, ClusterResList)
  #Order the dataframe by number of clusters
  Sankey_Full_Clust_df = Sankey_Full_Clust_df[order(grepl("_cluster", names(Sankey_Full_Clust_df)), order(gtools::mixedorder(names(Sankey_Full_Clust_df))))]
  #create a list of colors for Sankey plot
  Sankey_Colors <- RColorBrewer::brewer.pal(9, "Set1")
  Sankey_Link_List = list()
  #Create source-target dataframes for each set of cluster numbers
  for(i in 1:(ncol(Sankey_Full_Clust_df)-1)){
    #save the column name of both columns
    clustname_i = names(Sankey_Full_Clust_df)[i]
    clustname_j = names(Sankey_Full_Clust_df)[i+1]
    #save the number of clusters assessed in method
    clustnum_i = strsplit(clustname_i, split = "_")[[1]][2]
    clustnum_j = strsplit(clustname_j, split = "_")[[1]][2]
    #create a name so each dataframe can be identified
    Sankey_name = paste0(clustnum_i, "_", clustnum_j)
    Sankey_Single_Link = Sankey_Full_Clust_df |>
      dplyr::select(c(dplyr::all_of(i), dplyr::all_of(i + 1))) |> #select only the source and target column
      dplyr::mutate("{clustname_i}" := paste0(Sankey_Full_Clust_df[[i]], " (cl", clustnum_i, ")")) |> #add a name indicating the cluster number for source column
      dplyr::mutate("{clustname_j}" := paste0(Sankey_Full_Clust_df[[i+1]], " (cl", clustnum_j, ")")) |> #add a name indicating the cluster number for source column
      dplyr::select(source = dplyr::all_of(clustname_i), target = dplyr::all_of(clustname_j)) |> #rename the columns to source and target
      dplyr::mutate(value = 1L) #add a value column
    Sankey_Single_Link = cbind(Sankey_Single_Link, Metadata = as.factor(MetadataColumn), MetadataColor = as.factor(MetadataColumn)) #bind the metadata and color columns
    levels(Sankey_Single_Link$MetadataColor) = Sankey_Colors[1:nlevels(Sankey_Single_Link$MetadataColor)] #use the Sankey_Colors list to select colors for the data based on metadata levels
    Sankey_Single_Link = tibble::rownames_to_column(Sankey_Single_Link, "Sample_ID") #create a Sample_ID column
    Sankey_Link_List[[Sankey_name]] = Sankey_Single_Link #Save in a list
  }
  Sankey_Links = dplyr::bind_rows(Sankey_Link_List) #bind rows from all dataframes in list together
  Sankey_Nodes = data.frame(name = unique(c(Sankey_Links$source, Sankey_Links$target))) #create the node dataframe by identifying unique sources and targets
  Sankey_Nodes$color = as.factor(c("grey")) #create the node dataframe by identifying unique sources and targets
  Sankey_Links$source_id <- match(Sankey_Links$source, Sankey_Nodes$name) - 1
  Sankey_Links$target_id <- match(Sankey_Links$target, Sankey_Nodes$name) - 1
  Sankey_List = list(Links = Sankey_Links, Nodes = Sankey_Nodes) #create a list of the Links and Nodes dataframe
  return(Sankey_List)
}

#Sankey plot from Sankey data manipulated list
#' Sankey plot for sample movement throughout clusters
#' @description
#' Using the list of dataframes created in either
#' `Data_Manipulation_Sankey_Clusters()` or `Data_Manipulation_Sankey_Methods()`,
#' create a Sankey plot displaying the sample movement through different cluster
#' assignments. An interactive Sankey plot with samples colored by a metadata
#' feature.
#' @param SankeyList List of Sankey plot input, as created in either
#' `Data_Manipulation_Sankey_Clusters()` or `Data_Manipulation_Sankey_Methods()`.
#' @param Title String indicating the title for the plot. Default is
#' `"Sankey Plot"`.
#' @returns Sankey plot displaying samples through different cluster assignments,
#' with samples colored for a metadata feature, and a legend of this metadata
#' feature.
#' @keywords internal
Sankey_Plot = function(SankeyList,
                       Title = "Sankey Plot"){
  #create legend color vector
  Sankey_Legend_df = SankeyList$Links[, c("Metadata", "MetadataColor")]
  Sankey_Legend_df = Sankey_Legend_df[!duplicated(Sankey_Legend_df),]

  rownames(Sankey_Legend_df) = Sankey_Legend_df$Metadata
  Sankey_Legend_df$Metadata = NULL
  Sankey_Color = Sankey_Legend_df[, "MetadataColor"]
  names(Sankey_Color) = rownames(Sankey_Legend_df)

  #create legend scatterplot
  Sankey_Legend = plotly::plot_ly(type = "scatter",
                                  mode = "markers")
  for(Sample in names(Sankey_Color)){
    Sankey_Legend = plotly::add_trace(
      p = Sankey_Legend,
      x = 0,
      y = 0,
      name = Sample,
      mode = "lines",
      line = list(color = as.vector(Sankey_Color[Sample]), width = 10),
      showlegend = TRUE
    )
  }
  Sankey_Legend = plotly::layout(
    p = Sankey_Legend,
    xaxis = list(visible = FALSE, range = c(1,2)),
    yaxis = list(visible = FALSE, range = c(1,2)),
    showlegend = TRUE
  )

  Plot = plotly::plot_ly(type = "sankey",
                         domain = list(x = c(0, 0.95), y = c(0, 1)),
                         orientation = "h",
                         node = list(label = SankeyList$Nodes$name,
                                     color = SankeyList$Nodes$color),
                         link = list(source = SankeyList$Links$source_id,
                                     target = SankeyList$Links$target_id,
                                     value = SankeyList$Links$value,
                                     label = SankeyList$Links$Sample_ID,
                                     color = SankeyList$Links$MetadataColor))
  Plot = plotly::layout(p = Plot,
                        title = list(text = Title))

  Sankey_Subplot = plotly::subplot(
    Plot,
    Sankey_Legend,
    nrows = 1,
    margin = 0.05,
    widths = c(0.8, 0.2)
  )

  return(Sankey_Subplot)
}

#Data manipulation for the Sankey plot comparing different methods results
#' Create dataframe to be used for Sankey plot comparing multi-omics integration
#' methods
#' @description
#' Combine the cluster results for multiple multi-omics methods performed for
#' a single cluster number in a single dataframe. For each sample with each
#' different cluster assignment, the sample is given a source and target value,
#' which can be used to display the flow in the Sankey plot. Each sample also
#' gets a color assigned based on a user-defined metadata feature.
#' @param ClusterResData List of cluster assignment dataframes.
#' @param MetadataColumn String of the sample metadata feature for which the
#' samples should be colored. The metadata feature needs to be available
#' in the metadata files for all omics types.
#' @returns List with results. $Links, dataframe with samples in rows, and
#' columns indicating source cluster and target cluster, as well as metadata
#' color. $Nodes, dataframe with source and target clusters.
#' @keywords internal
Data_Manipulation_Sankey_Methods = function(ClusterResData,
                                            MetadataColumn){
  #create list of colors for Sankey plot
  Sankey_Colors_Short = RColorBrewer::brewer.pal(9, "Set1")
  Sankey_Colors_Long = Polychrome::glasbey.colors(n = 32)
  Sankey_Greys = RColorBrewer::brewer.pal(9, "Greys")
  Sankey_Link_List_Method = list()
  #create source-target dataframes for each set of methods
  for(i in 1:(ncol(ClusterResData)-1)){
    #save the column name of both columns
    Sankey_clustname_i = names(ClusterResData)[i]
    Sankey_clustname_j = names(ClusterResData)[i+1]
    #save abbreviation of name for plot
    if(grepl("MoCluster", Sankey_clustname_i)){
      Sankey_methabbr_i = "MoCl"
    } else if(grepl("MCIA", Sankey_clustname_i)){
      Sankey_methabbr_i = "MCIA"
    } else if(grepl("jNMF", Sankey_clustname_i)){
      Sankey_methabbr_i = "jNMF"
    } else if(grepl("iNMF", Sankey_clustname_i)){
      Sankey_methabbr_i = "iNMF"
    } else if(grepl("LRAcluster", Sankey_clustname_i)){
      Sankey_methabbr_i = "LRA"
    } else if(grepl("COCA", Sankey_clustname_i)){
      Sankey_methabbr_i = "COCA"
    } else if(grepl("MOFA", Sankey_clustname_i)){
      Sankey_methabbr_i = "MOFA"
    } else if(grepl("GAUDI", Sankey_clustname_i)){
      if(grepl("HDBSCAN", Sankey_clustname_i)){
        Sankey_methabbr_i = "GAUDI_HDB"
      } else if(grepl("DBSCAN", Sankey_clustname_i)){
        Sankey_methabbr_i = "GAUDI_DB"
      } else if(grepl("Kmeans", Sankey_clustname_i)){
        Sankey_methabbr_i = "GAUDI_K"
      }
    } else if(grepl("SNF", Sankey_clustname_i)){
      if(grepl("euclidean squared", Sankey_clustname_i)){
        Sankey_methabbr_i = "SNF_EuSq"
      } else if(grepl("euclidean", Sankey_clustname_i)){
        Sankey_methabbr_i = "SNF_Eu"
      } else if(grepl("manhattan", Sankey_clustname_i)){
        Sankey_methabbr_i = "SNF_Man"
      } else if(grepl("minkowski 0.25", Sankey_clustname_i)){
        Sankey_methabbr_i = "SNF_Min_025"
      } else if(grepl("minkowski 0.5", Sankey_clustname_i)){
        Sankey_methabbr_i = "SNF_Min_05"
      } else if(grepl("minkowski 3", Sankey_clustname_i)){
        Sankey_methabbr_i = "SNF_Min_3"
      } else if(grepl("minkowski 4", Sankey_clustname_i)){
        Sankey_methabbr_i = "SNF_Min_4"
      }
    } else if(grepl("EnsembleCHC", Sankey_clustname_i)){
      Sankey_methabbr_i = "CHC"
    } else if(grepl("EnsembleCCC", Sankey_clustname_i)){
      Sankey_methabbr_i = "CCC"
    } else if(grepl("EnsembleMDS_HC", Sankey_clustname_i)){
      Sankey_methabbr_i = "MDS_HC"
    } else if(grepl("EnsembleMDS_km", Sankey_clustname_i)){
      Sankey_methabbr_i = "MDS_km"
    } else if(grepl("EnsembleMPA", Sankey_clustname_i)){
      Sankey_methabbr_i = "MPA"
    } else if(grepl("EnsembleCCA", Sankey_clustname_i)){
      Sankey_methabbr_i = "CCA"
    } else if(grepl("EnsembleCAT", Sankey_clustname_i)){
      Sankey_methabbr_i = "CAT"
    } else if(grepl("EnsembleUMAPCCA", Sankey_clustname_i)){
      Sankey_methabbr_i = "CCA_UMAP"
    } else if(grepl("EnsembleUMAPCAT", Sankey_clustname_i)){
      Sankey_methabbr_i = "CAT_UMAP"
    }
    if(grepl("MoCluster", Sankey_clustname_j)){
      Sankey_methabbr_j = "MoCl"
    } else if(grepl("MCIA", Sankey_clustname_j)){
      Sankey_methabbr_j = "MCIA"
    } else if(grepl("jNMF", Sankey_clustname_j)){
      Sankey_methabbr_j = "jNMF"
    } else if(grepl("iNMF", Sankey_clustname_j)){
      Sankey_methabbr_j = "iNMF"
    } else if(grepl("LRAcluster", Sankey_clustname_j)){
      Sankey_methabbr_j = "LRA"
    } else if(grepl("COCA", Sankey_clustname_j)){
      Sankey_methabbr_j = "COCA"
    } else if(grepl("MOFA", Sankey_clustname_j)){
      Sankey_methabbr_j = "MOFA"
    } else if(grepl("GAUDI", Sankey_clustname_j)){
      if(grepl("HDBSCAN", Sankey_clustname_j)){
        Sankey_methabbr_j = "GAUDI_HDB"
      } else if(grepl("DBSCAN", Sankey_clustname_j)){
        Sankey_methabbr_j = "GAUDI_DB"
      } else if(grepl("Kmeans", Sankey_clustname_j)){
        Sankey_methabbr_j = "GAUDI_K"
      }
    } else if(grepl("SNF", Sankey_clustname_j)){
      if(grepl("euclidean squared", Sankey_clustname_j)){
        Sankey_methabbr_j = "SNF_EuSq"
      } else if(grepl("euclidean", Sankey_clustname_j)){
        Sankey_methabbr_j = "SNF_Eu"
      } else if(grepl("manhattan", Sankey_clustname_j)){
        Sankey_methabbr_j = "SNF_Man"
      } else if(grepl("minkowski 0.25", Sankey_clustname_j)){
        Sankey_methabbr_j = "SNF_Min_025"
      } else if(grepl("minkowski 0.5", Sankey_clustname_j)){
        Sankey_methabbr_j = "SNF_Min_05"
      } else if(grepl("minkowski 3", Sankey_clustname_j)){
        Sankey_methabbr_j = "SNF_Min_3"
      } else if(grepl("minkowski 4", Sankey_clustname_j)){
        Sankey_methabbr_j = "SNF_Min_4"
      }
    } else if(grepl("EnsembleCHC", Sankey_clustname_j)){
      Sankey_methabbr_j = "CHC"
    } else if(grepl("EnsembleCCC", Sankey_clustname_j)){
      Sankey_methabbr_j = "CCC"
    } else if(grepl("EnsembleMDS_HC", Sankey_clustname_j)){
      Sankey_methabbr_j = "MDS_HC"
    } else if(grepl("EnsembleMDS_km", Sankey_clustname_j)){
      Sankey_methabbr_j = "MDS_km"
    } else if(grepl("EnsembleMPA", Sankey_clustname_j)){
      Sankey_methabbr_j = "MPA"
    } else if(grepl("EnsembleCCA", Sankey_clustname_j)){
      Sankey_methabbr_j = "CCA"
    } else if(grepl("EnsembleCAT", Sankey_clustname_j)){
      Sankey_methabbr_j = "CAT"
    } else if(grepl("EnsembleUMAPCCA", Sankey_clustname_j)){
      Sankey_methabbr_j = "CCA_UMAP"
    } else if(grepl("EnsembleUMAPCAT", Sankey_clustname_j)){
      Sankey_methabbr_j = "CAT_UMAP"
    }
    #create a name to identify each source-target dataframe by methods
    Sankey_name = paste0(Sankey_methabbr_i, "_", Sankey_methabbr_j)
    Sankey_Single_Link_Method = ClusterResData |>
      dplyr::select(c(dplyr::all_of(i), dplyr::all_of(i + 1))) |> #select only the source and target column
      dplyr::mutate("{Sankey_clustname_i}" := paste0(ClusterResData[[i]], " (", Sankey_methabbr_i, ")")) |> #add a name indicating the cluster number for source column
      dplyr::mutate("{Sankey_clustname_j}" := paste0(ClusterResData[[i+1]], " (", Sankey_methabbr_j, ")")) |> #add a name indicating the cluster number for source column
      dplyr::select(source = dplyr::all_of(Sankey_clustname_i), target = dplyr::all_of(Sankey_clustname_j)) |> #rename the columns to source and target
      dplyr::mutate(value = 1L) #add a value column
    Sankey_Single_Link_Method = cbind(Sankey_Single_Link_Method, Metadata = as.factor(MetadataColumn), MetadataColor = as.factor(MetadataColumn)) #bind the metadata and color columns
    if(nlevels(Sankey_Single_Link_Method$MetadataColor) <= 10){
      Sankey_Colors = Sankey_Colors_Short
    }else{
      Sankey_Colors = Sankey_Colors_Long
    }
    levels(Sankey_Single_Link_Method$MetadataColor) = Sankey_Colors[1:nlevels(Sankey_Single_Link_Method$MetadataColor)] #use the Sankey_Colors list to select colors for the data based on metadata levels
    Sankey_Single_Link_Method = tibble::rownames_to_column(Sankey_Single_Link_Method, "Sample_ID") #create a Sample_ID column
    Sankey_Link_List_Method[[Sankey_name]] = Sankey_Single_Link_Method #Save in a list
  }
  Sankey_Links_Method = dplyr::bind_rows(Sankey_Link_List_Method) #bind rows from all dataframes in list together
  Sankey_Nodes_Method = data.frame(name = unique(c(Sankey_Links_Method$source, Sankey_Links_Method$target))) #create the node dataframe by identifying unique sources and targets
  Sankey_Nodes_Method$color = as.factor(c("grey"))
  Sankey_Links_Method$source_id <- match(Sankey_Links_Method$source, Sankey_Nodes_Method$name) - 1
  Sankey_Links_Method$target_id <- match(Sankey_Links_Method$target, Sankey_Nodes_Method$name) - 1
  Sankey_List_Method = list(Links = Sankey_Links_Method, Nodes = Sankey_Nodes_Method) #create a list of the Links and Nodes dataframe
  return(Sankey_List_Method)
}

#create top feature weights plot per dataframe
#' Plot features with heightest absolute weigths per multi-omics integration
#' method
#' @description
#' Plot the features with the highest absolute weights in each factor
#' (component), per omics type and multi-omics integration method.
#' @param FeatureRes Dataframe with feature weights per factor (component).
#' @param NumberFeatures Number of features to plot for each factor. Default is
#' 10.
#' @param Scale Whether to scale the data from the different omics types, to
#' make the weights more comparable between omics types. Default is `TRUE`.
#' @returns Plot with the top n absolute features by weight, for each factor
#' (component).
#' @keywords internal
Top_Feature_Weights_All_Factors = function(FeatureRes,
                                           NumberFeatures = 10,
                                           Scale = TRUE){
  FeatureResDF = tibble::rownames_to_column(FeatureRes, "Features")
  LongerDF = FeatureResDF |> tidyr::pivot_longer(!Features, names_to = "Factor", values_to = "Weight")
  LongerDF$Factor = as.factor(LongerDF$Factor)
  if(Scale == TRUE){
    LongerDF = LongerDF |> dplyr::mutate(Weight = Weight/max(abs(Weight)))
  }
  LongerDFTop = LongerDF |> dplyr::group_by(Factor) |> dplyr::slice_max(abs(Weight), n = NumberFeatures) |> dplyr::ungroup()
  LongerDFTop$AbsWeight = abs(LongerDFTop$Weight)
  LongerDFTop$Sign = ifelse(LongerDFTop$Weight > 0, "+", "-")
  if(any(duplicated(LongerDFTop$Features))){
    LongerDFTop$Features = paste0(LongerDFTop$Features, "_F", substring(LongerDFTop$Factor, 7))
  }
  LongerDFTop$Features = factor(LongerDFTop$Features, levels = rev(unique(LongerDFTop$Features)))

  Plot = ggplot2::ggplot(LongerDFTop, ggplot2::aes(x = Features, y = AbsWeight)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_segment(ggplot2::aes(xend = Features), linewidth = 0.75, yend = 0) +
    ggplot2::coord_flip() +
    ggplot2::labs(y = "Weight") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_text(color = "black"),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(color = "black"),
                   axis.text.y = ggplot2::element_text(size = ggplot2::rel(1.1), hjust = 1, color = "black"),
                   axis.ticks.y = ggplot2::element_line(),
                   strip.text = ggplot2::element_text(size = ggplot2::rel(1.2)),
                   panel.background = ggplot2::element_blank(),
                   panel.spacing = ggplot2::unit(1, "lines"),
                   panel.grid.major.y = ggplot2::element_blank()) +
    ggplot2::facet_wrap(~Factor, ncol = 3, scales = "free")

  Plot = Plot +
    ggplot2::ylim(0, max(LongerDFTop$AbsWeight) + 0.1) +
    ggplot2::geom_text(label = LongerDFTop$Sign, y = max(LongerDFTop$AbsWeight) + 0.1, size = 5)

  return(Plot)
}

#dendrogram metadata coloring
#' Create list of dataframes with sample metadata information and color for
#' colored bars below dendrograms
#' @description
#' Create list of dataframes with sample metadata information and colors. These
#' can be used for create colored bars and colored labels on sample dendrograms.
#' @param MetadataFeatures Names of the sample metadata features to be
#' included in colored bars. A string vector with the names of the metadata
#' features as in the metadata files should be provided.
#' @param Metadata Dataframe with sample metadata.
#' @param Dendrogram Sample dendrogram.
#' @param LabelMetadata Name of the metadata feature for which the sample
#' labels should be colored. Must be a sting with the name of the metadata
#' feature as in the metadata files.
#' @returns List of dataframes. $Vector, dataframe with colored levels of
#' metadata features for the colored bars. $Label, dataframe with labels colored
#' by metadata features, ordered by the dendrogram order.
#' @keywords internal
Dendrogram_Sample_Meta = function(MetadataFeatures,
                                  Metadata,
                                  Dendrogram,
                                  LabelMetadata){
  Sample_Metadata = list()
  Vector = list()
  Meta_Short = list(Meta1 = RColorBrewer::brewer.pal(9, "Set1"),
                    Meta2 = RColorBrewer::brewer.pal(8, "Set2"),
                    Meta3 = RColorBrewer::brewer.pal(8, "Dark2"),
                    Meta4 = RColorBrewer::brewer.pal(8, "Accent"),
                    Meta5 = RColorBrewer::brewer.pal(12, "Set3"),
                    Meta6 = RColorBrewer::brewer.pal(9, "Pastel1"))
  Meta_Long = list(Meta1 = Polychrome::glasbey.colors(n = 32),
                   Meta2 = Polychrome::palette36.colors(n = 36),
                   Meta3 = Polychrome::dark.colors(n = 24),
                   Meta4 = Polychrome::alphabet.colors(n = 26),
                   Meta5 = Polychrome::sky.colors(n = 24),
                   Meta6 = Polychrome::green.armytage.colors(n = 26))

  for(i in 1:length(MetadataFeatures)){
    Feature = MetadataFeatures[i]
    Meta_Factor = as.factor(Metadata[, Feature])
    if(nlevels(Meta_Factor) <= 10){
      Meta_Color = Meta_Short[[i]]
    }else{
      Meta_Color = Meta_Long[[i]]
    }
    Meta_Colored = Meta_Color[Meta_Factor]
    Meta_Vector = as.vector(Meta_Colored)
    Vector[[Feature]] = Meta_Vector
    if(Feature == LabelMetadata){
      Meta_Order = Meta_Colored[stats::order.dendrogram(Dendrogram)]
      Sample_Metadata$Label = Meta_Order
    }
  }
  Vector_df = as.data.frame(do.call(cbind, Vector))
  Vector_df = Vector_df |> dplyr::relocate(dplyr::all_of(LabelMetadata), .after = last_col())
  Sample_Metadata$Vector = Vector_df
  return(Sample_Metadata)
}

#dendrogram label coloring
#' Create a list of dataframes with label coloring and legend information for a
#' dendrogram
#' @description
#' Create a list of dataframes with dendrogram labels colored and legend
#' information for the label colors. To be used with a dendrogram.
#' @param Feature Whether feature or factor dendrogram is being used. If `TRUE`,
#' the input dendrogram is a feature dendrogram. If `FALSE`, the input
#' dendrogram is a factor dendrogram.
#' @param HClust Hierarchical clustering trees.
#' @param Dendrogram Dendrogram.
#' @param LabelSize Whether to alternate label size in feature dendrogram, to
#' distinguish the omics.
#' @returns List of dataframes. $Colors, dataframe with dendrogram labels and
#' the colors to be used. $Legend_Levels, dataframe with dendrogram label levels.
#' $Legend_Colors, dataframe with dendrogram label level colors.
#' @keywords internal
Dendrogram_Label_Groups = function(Feature,
                                   HClust,
                                   Dendrogram,
                                   LabelSize){
  FeatureNames = HClust$labels
  if(Feature == T){
    Split_DF = data.frame(do.call(rbind, strsplit(FeatureNames, "_\\s*(?=[^_]+$)", perl = T)))
    Split_DF$X3 = FeatureNames
    Groups = as.factor(Split_DF$X2)
  }else if(Feature == F){
    Split_DF = data.frame(do.call(rbind, strsplit(FeatureNames, "_")))
    Split_DF$X3 = FeatureNames
    Groups = as.factor(Split_DF$X1)
  }
  Colors = RColorBrewer::brewer.pal(length(Groups), name = "Set1")
  Colored = Colors[Groups]
  Ordered = Colored[stats::order.dendrogram(Dendrogram)]

  Fixed_Labels = list()
  Fixed_Labels$Colors = Ordered
  Fixed_Labels$Legend_Levels = Groups
  Fixed_Labels$Legend_Colors = Colors

  if(Feature == T & LabelSize == T){
    Levels = levels(Groups)
    n = round(length(Levels)/2, 0)
    Group1 = sample(Levels, n)
    Group2 = Levels[!(Levels %in% Group1)]
    Split_DF = Split_DF |>
      dplyr::rowwise() |>
      dplyr::mutate(Labels = dplyr::case_when(X2 %in% Group1 ~ X3,
                                              X2 %in% Group2 ~ X1)) |>
      dplyr::ungroup()
    Sized_Labels = Split_DF$Labels
    Order_Labels = Sized_Labels[stats::order.dendrogram(Dendrogram)]
    Fixed_Labels$Sized = Order_Labels
  }
  return(Fixed_Labels)
}

#dendrogram omics type colored bar
#' Create list of dataframes with data for colored bars and legend for
#' dendrogram
#' @description
#' Create a list of dataframes with omics information for the colored bars and
#' legend data for feature dendrograms.
#' @param HClust Hierachical clustering tree.
#' @param Dendrogram Dendrogram.
#' @returns A list. $Colors, dataframe with the color to be used per feature
#' label. $Legend_Levels, dataframe with dendrogram omics levels. $Legend_Colors,
#' dataframe with dendrogram omics level colors. $Vector, vector with color to
#' be used per feature, ordered by dendrogram order.
#' @keywords internal
Dendrogram_Omics_Bar = function(HClust,
                                Dendrogram){
  FeatureNames = HClust$labels
  Split_DF = data.frame(do.call(rbind, strsplit(FeatureNames, "_\\s*(?=[^_]+$)", perl = T)))
  Omics = as.factor(Split_DF$X2)

  Colors = RColorBrewer::brewer.pal(length(Omics), name = "Set1")
  Colored = Colors[Omics]
  Ordered = Colored[stats::order.dendrogram(Dendrogram)]
  Vector = as.vector(Colored)

  Omics_Bar = list()
  Omics_Bar$Colors = Ordered
  Omics_Bar$Legend_Levels = Omics
  Omics_Bar$Legend_Colors = Colors
  Omics_Bar$Vector = Vector

  return(Omics_Bar)
}

#----------------EINS R6 Merge----------------
#' Merge data from one EINS R6 object into another EINS R6 object
#' @description
#' Merge the data from one EINS R6 object with data from another EINS R6 object.
#' Can be used if different preprocessing is wanted for separate multi-omics
#' integration methods, or when performing computations in parallel. Preprocessed
#' omics, single omics clustering, multi-omics integration or ensemble integration
#' results can all be merged. If data is present in both EINS R6 objects, data in
#' the reference EINS R6 object will be overwritten with data from the
#' supplemental EINS R6 object.
#' @param Ref_EINS  Reference EINS R6 object, from which all data will be taken.
#' @param Supp_EINS Supplemental EINS R6 object, from which only specified data
#' will be taken, overwriting data in the reference EINS R6 object.
#' @param Omics List of omics names. For each omics type provided in the list,
#' the preprocessed matrices will be stored in the new EINS R6 object. If the
#' reference EINS R6 contains preprocessed omics matrices for the specified omics,
#' these will be overwritten.
#' @param Single_Omics List of omics names. For each omics type provided in the
#' list, the single omics clustering results will be stored in the new EINS R6
#' object. If the reference EINS R6 contains single omics clustering results for
#' the specified omics, these will be overwritten.
#' @param Multi_Omics List of multi-omics integration methods. For each multi-omics
#' integration method provided in the list, the integration results will be
#' stored in the new EINS R6 object. If the reference EINS R6 contains multi-omics
#' integration results for the specified methods, these will be overwritten.
#' @param Ensemble List of ensemble integration methods. For each ensemble
#' integration method provided in the list, the ensemble integration results will
#' be stored in the new EINS R6 object. If the reference EINS R6 contains ensemble
#' integration results for the specified methods, these will be overwritten.
#' @returns An EINS R6 object. Contains the data from the reference EINS R6 object,
#' except for the supplied omics and methods, which are taken from the supplemental
#' EINS R6 object.
EINS_Merge = function(Ref_EINS = NULL,
                      Supp_EINS = NULL,
                      Omics = NA,
                      Single_Omics = NA,
                      Multi_Omics = NA,
                      Ensemble = NA
){
  EINS_3 = Ensemble_Integration_Suite$new()
  EINS_3$Omics = Ref_EINS$Omics
  EINS_3$Data_Checks = Ref_EINS$Data_Checks
  EINS_3$Preprocessed_Omics = Ref_EINS$Preprocessed_Omics
  EINS_3$Single_Omics = Ref_EINS$Single_Omics
  EINS_3$Multi_Omics = Ref_EINS$Multi_Omics
  EINS_3$Ensemble = Ref_EINS$Ensemble

  for(Omic in Omics){
    if(is.na(Omics)){
      next
    }else{
      EINS_3$Preprocessed_Omics$Data[[Omic]] = Supp_EINS$Preprocessed_Omics$Data[[Omic]]
      EINS_3$Preprocessed_Omics$Matched_Data[[Omic]] = Supp_EINS$Preprocessed_Omics$Matched_Data[[Omic]]
      EINS_3$Preprocessed_Omics$Non_Neg_Data[[Omic]] = Supp_EINS$Preprocessed_Omics$Non_Neg_Data[[Omic]]
      EINS_3$Preprocessed_Omics$Correlated_Features[[Omic]] = Supp_EINS$Preprocessed_Omics$Correlated_Features[[Omic]]
    }
  }

  for(Omic in Single_Omics){
    if(is.na(Omic)){
      next
    }else{
      EINS_3$Single_Omics$Fit[[Omic]] = Supp_EINS$Single_Omics$Fit[[Omic]]
      EINS_3$Single_Omics$ClusterRes[[Omic]] = Supp_EINS$Single_Omics$ClusterRes[[Omic]]
      EINS_3$Single_Omics$HClustTree[[Omic]] = Supp_EINS$Single_Omics$HClustTree[[Omic]]
      EINS_3$Single_Omics$DistanceMatrix[[Omic]] = Supp_EINS$Single_Omics$DistanceMatrix[[Omic]]
    }
  }

  for(MO in Multi_Omics){
    if(is.na(MO)){
      next
    }else{
      EINS_3$Multi_Omics$Fit[[MO]] = Supp_EINS$Multi_Omics$Fit[[MO]]
      EINS_3$Multi_Omics$ClusterRes[[MO]] = Supp_EINS$Multi_Omics$ClusterRes[[MO]]
      EINS_3$Multi_Omics$FeatureRes[[MO]] = Supp_EINS$Multi_Omics$FeatureRes[[MO]]
      EINS_3$Multi_Omics$FactorRes[[MO]] = Supp_EINS$Multi_Omics$FactorRes[[MO]]
      EINS_3$Multi_Omics$VarianceExp[[MO]] = Supp_EINS$Multi_Omics$VarianceExp[[MO]]
      EINS_3$Multi_Omics$CoordData[[MO]] = Supp_EINS$Multi_Omics$CoordData[[MO]]
      EINS_3$Multi_Omics$Feature_HClust$DistanceMatrix[[MO]] = Supp_EINS$Multi_Omics$Feature_HClust$DistanceMatrix[[MO]]
      EINS_3$Multi_Omics$Feature_HClust$HClustTree[[MO]] = Supp_EINS$Multi_Omics$Feature_HClust$HClustTree[[MO]]
    }
  }

  for(ENS in Ensemble){
    if(is.na(ENS)){
      next
    }else{
      EINS_3$Ensemble$Samples[[ENS]] = Supp_EINS$Ensemble$Samples[[ENS]]
      EINS_3$Ensemble[[ENS]] = Supp_EINS$Ensemble[[ENS]]
    }
  }
  return(EINS_3)
}

#----------------END----------------


