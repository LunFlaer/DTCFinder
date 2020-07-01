#' FilterData Function loads in data from 10X
#'
#' This function allows you to input 10X single-cell RNAseq profile
#' (loading sparse data matrices provided by 10X genomics),
#' and remove some unwanted cells.
#'
#' @param PathToData Folder with barcodes.tsv, features.tsv (or genes.tsv), matrix.mtx. Default is the current folder
#' @param FeatureColumn choose the column of genes.tsv or features.tsv for feature identifiers; default is 2
#' @param UniqueFeatures Logical for unique feature identifiers (default TRUE)
#' @param FeatureSuffixTrim Logicalfor removal of trailing "-1" in cell barcodes (default FALSE)
#' @param MinFeatures Threshold to remove those cell barcodes with found features fewer than it (default 200)
#' @param MinUMIs threshold to remove those cell barcodes with found UMIs fewer than it (default 500)
#'
#' @return a sparse matrix containing the expression data will be returned.
#'
#' @importFrom Matrix readMM
#'
#' @export
#' @examples
#' # For output from CellRanger < 3.0
#' FilterData("outs/filtered_gene_bc_matrices_mex/GRCh38/")
#' # For output from CellRanger >= 3.0
#' FilterData("outs/filtered_feature_bc_matrix/")
#'

FilterData <- function(
  PathToData = ".",
  FeatureColumn = 2,
  UniqueFeatures = TRUE,
  FeatureSuffixTrim = FALSE,
  MinFeatures = 200,
  MinUMIs = 500
){
  if (!dir.exists(paths = PathToData)) {
      stop(paste0("Error: Folder '", PathToData, "' does not exist"))
  }
  features.file <- file.path(PathToData, 'genes.tsv')
  barcodes.file <- file.path(PathToData, 'barcodes.tsv')
  matrix.file <- file.path(PathToData, 'matrix.mtx')
  if (!file.exists(features.file)) {
    features.file <- file.path(PathToData, 'features.tsv.gz')
    barcodes.file <- paste0(barcodes.file, ".gz")
    matrix.file <- paste0(matrix.file, ".gz")
  }
  if (!file.exists(barcodes.file)) {
    stop(paste0("Error: Barcode file '", basename(path = barcodes.file), "' missing."))
  }
  if (!file.exists(features.file) ) {
    stop(paste0("Error: Features file '", basename(path = features.file), "' missing"))
  }
  if (!file.exists(matrix.file)) {
    stop(paste0("Error: Expression matrix file '", basename(path = matrix.file), "' missing."))
  }
  if (!require("Matrix", quietly = TRUE)){
    install.packages("Matrix")
  }
  data <- Matrix::readMM(file = matrix.file)	# MatrixMarket format
  barcodes <- readLines(barcodes.file)
  if (all(grepl(pattern = "\\-1$", x = barcodes)) & FeatureSuffixTrim) {
    barcodes <- unlist(strsplit(x = barcodes, split = "-", fixed=TRUE))[which(nchar(unlist(strsplit(x = barcodes, split = "-", fixed=TRUE)))>2)]
  }
  colnames(data) <- barcodes
  features <- read.delim(file = features.file, header = FALSE, stringsAsFactors = FALSE)
  if (any(is.na(x = features[, FeatureColumn]))) {
      warning(
        'There are some NA features. Replacing NA features with IDs from another column now',
        call. = FALSE,
        immediate. = TRUE
      )
      na.features <- which(is.na(features[, FeatureColumn]))
      replacement.column <- ifelse(test = FeatureColumn == 2, yes = 1, no = 2)
      features[na.features, FeatureColumn] <- features[na.features, replacement.column]
  }
  if (UniqueFeatures) {
      if (FeatureColumn > ncol(features)) {
        stop(paste0("feature column ", FeatureColumn,
                    " larger than the column number of feature.tsv.gz (or genes.tsv) ", ncol(features), ".",
                    " Try a smaller FeatureColumn argument ( <= ", ncol(features), " )."))
      }
      rownames(data) <- make.unique(features[, FeatureColumn])
  }
  if (ncol(features) > 2) {
    data = data[features[,3] == "Gene Expression",]
  }
  nFeatures = apply(data, 1, function(x) {sum(x>0, na.rm=TRUE)})
  nUMIs = apply(data, 2, sum)  # colSums(data_mat)
  barcodes_selected = which(nFeatures > MinFeatures & nUMIs > MinUMIs)
  return(data[barcodes_selected,])
}

