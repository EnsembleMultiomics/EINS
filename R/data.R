#' NCI60 toy data.
#'
#' List of reduced NCI60 proteomics, transcriptomics and DNA methylation data,
#' with metadata. Omics data is available for 56 cancer cell lines. 300 proteins
#' were selected (150 top proteins by standard deviation, 150 randomly
#' selected). 400 methylated genes were selected (200 top genes by standard
#' deviation, 200 randomly selected). 500 transcripts were selected (250 top
#' transcripts by standard deviation, 250 randomly selected). Metadata is
#' available for all 56 cell lines.
#'
#' @format A list with four underlying dataframes:
#' \describe{
#'  \dataframe{Proteomics}{300 rows of proteins, 56 columns of NCI60 cancer
#'  cell lines},
#'  \dataframe{Methylation}{400 rows of methylated genes, 56 columns of NCI60
#'  cancer cell lines},
#'  \dataframe{Transcriptomics}{500 rows of gene transcripts, 56 columns of
#'  NCI60 cancer cell lines},
#'  \dataframe{Metadata}{56 rows of NCI60 cancer cell lines, 14 columns of
#'  metadata features}
#'  \describe{
#'   \item{tissue of origin}{tissue of origin for the tumour cell line
#'   (breast, central nervous system, colon, leukemia, melanoma, non-small cell
#'   lung, ovarian, renal)},
#'   \item{age}{patient age in years at time of tumour sampling (4 - 75)},
#'   \item{sex}{patient sex},
#'   \item{prior treatment}{treatment received prior to tumour sampling},
#'   \item{Epithelial}{epithelial nature of the tissue of origin},
#'   \item{histology}{histological features of tumour},
#'   \item{source}{source of sampled tumour},
#'   \item{ploidy}{ploidy information of the tumour},
#'   \item{p53}{p53 mutation status},
#'   \item{mdr}{MDR function},
#'   \item{doubling time}{time needed for doubling of tumour},
#'   \item{Institute}{institute where tumour was sampled},
#'   \item{Contributor}{Contributor at institute},
#'   \item{Reference}{reference for original description of sample}
#'  }
#'}
"EINS_NCI60_Toy"
