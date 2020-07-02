#' Load gene-barcode matrix from a specific H5 file
#'
#' @param Path to the gene-barcode matrix H5 file (NOT the pipestance path)
#' @param genome The desired genome (e.g., 'hg19' or 'mm10')
#' @return A GeneBCMatrix sparse matrix where rows are genes and columns are cell-barcodes
#' @export
#' @import Matrix
#' @import rhdf5
#' @examples
#' \dontrun{
#' # Load hg19 from a filtered matrix h5
#' gene_bc_matrix <- get_matrix_from_h5("/home/user/cellranger_output/outs/filtered_gene_bc_matrices_h5.h5", "hg19")
#' }
GetMatrixFromH5 <- function(filename, genome=NULL) {
  if(!file.exists(filename)) {
    stop(sprintf("Could not find matrix H5 file: \n\t %s\n", filename))
  }
  
  genomes <- attributes(h5dump(filename, load=FALSE))$names
  
  if (is.null(genome)) {
    if (length(genomes) > 1) {
      stop(sprintf("Multiple genomes found; please specify one. \n Genomes present: %s",paste(genomes, collapse=", ")))
    }
    genome <- genomes[1]
  }
  
  if (!(genome %in% genomes)) {
    stop(sprintf("Genome %s not found in file. \n Genomes present: %s",genome,paste(genomes, collapse=", ")))
  }
  
  dset <- h5read(filename, name = genome)
  sparse_mat <- sparseMatrix(i = dset$indices + 1, p = dset$indptr, x = as.numeric(dset$data), dims = dset$shape, giveCsparse = FALSE)
  genes <- data.frame(id = dset$genes, symbol = dset$gene_names, row.names = dset$genes)
  barcodes <- data.frame(barcode = dset$barcodes, row.names = dset$barcodes)
  gbm <- newGeneBCMatrix(sparse_mat, genes, barcodes)
  H5close()
  return (gbm)
}
