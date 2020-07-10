#' Load gene-barcode matrix from a specific H5 file
#'
#' @param Path to the gene-barcode matrix H5 file (NOT the pipestance path)
#' @param genome The desired genome (e.g., 'hg19' or 'mm10')
#' @return A GeneBCMatrix sparse matrix where rows are genes and columns are cell-barcodes
#' @export
#' @import Matrix
#' @import rhdf5
#' @import stringi
#' 
#' @examples
#' \dontrun{
#' # Load hg19 from a filtered matrix h5
#' gene_bc_matrix <- get_matrix_from_h5("/home/user/cellranger_output/outs/filtered_gene_bc_matrices_h5.h5", "hg19")
#' }
GetFileFromH5 <- function(filename, genome=NULL, filetype = "all") {
  if (!require("stringi", quietly = TRUE)){
    install.packages("stringi")
  }
  if (!require("rhdf5", quietly = TRUE)){
    install.packages("rhdf5")
  }
  
  
  if(!file.exists(filename)) {
    stop(sprintf("Could not find matrix H5 file: \n\t %s\n", filename))
  }
  filedir <- substr(filename, 0, stri_length(filename) - stri_length(basename(filename)))
  
  genomes <- h5dump(filename)$matrix$features$genome
  
  if (is.null(genome)) {
    if (length(unique(genomes)) > 1) {
      genome <- readline(prompt=sprintf("Multiple genomes found; please specify one. \n Genomes present: %s \n",paste(genomes, collapse=", ")))
    }
    else{
      genome <- genomes[1]
    }
  }
  
  if (!(genome %in% genomes)) {
    stop(sprintf("Genome %s not found in file. \n Genomes present: %s",genome,paste(genomes, collapse=", ")))
  }

  
  dset <- h5read(filename, name = "matrix")
  
  sparse_mat <- Matrix::sparseMatrix(i = dset$indices + 1, p = dset$indptr, x = as.numeric(dset$data), dims = dset$shape, giveCsparse = FALSE)
  
  
  
  if (filetype == "matrix") {
    if (!require("Matrix", quietly = TRUE)){
      install.packages("Matrix")
    }
    if (!require("cellrangerRkit", quietly = TRUE)){
      if (!require("BiocManager", quietly = TRUE)){
        install.packages("BiocManager")
      }
      devtools::install_github("hb-gitified/cellrangerRkit")
    }
    barcodes <- data.frame(barcode = dset$barcodes, row.names = NULL)
    genes <- data.frame(id = dset$features$id, symbol = dset$features$name, row.names = dset$features$id)
    gbm <- newGeneBCMatrix(sparse_mat, genes, barcodes)
    Matrix::writeMM(obj = sparse_mat, file=paste0(filedir, "matrix.mtx"))
  }
  else if (filetype == "gene") {
    genes <- data.frame(id = dset$features$id, symbol = dset$features$name, row.names = dset$features$id)
    write.table(genes, file=paste0(filedir, 'genes.tsv'), quote=FALSE, sep='\t', col.names = F, row.names = F)
  }
  else if (filetype == "barcode") {
    barcodes <- data.frame(barcode = dset$barcodes, row.names = NULL)
    write.table(barcodes, file=paste0(filedir, 'barcodes.tsv'), quote=FALSE, sep="\t", col.names = F, row.names = F)
  }
  else if (filetype == "all") {
    if (!require("Matrix", quietly = TRUE)){
      install.packages("Matrix")
    }
    if (!require("cellrangerRkit", quietly = TRUE)){
      if (!require("BiocManager", quietly = TRUE)){
        install.packages("BiocManager")
      }
      devtools::install_github("hb-gitified/cellrangerRkit")
    }
    
    barcodes <- data.frame(barcode = dset$barcodes, row.names = NULL)
    genes <- data.frame(id = dset$features$id, symbol = dset$features$name, row.names = dset$features$id)
    gbm <- newGeneBCMatrix(sparse_mat, genes, barcodes)
    
    Matrix::writeMM(obj = sparse_mat, file=paste0(filedir, "matrix.mtx"))
    write.table(genes, file=paste0(filedir, 'genes.tsv'), quote=FALSE, sep='\t', col.names = F, row.names = F)
    write.table(barcodes, file=paste0(filedir, 'barcodes.tsv'), quote=FALSE, sep="\t", col.names = F, row.names= F)
  }
  H5close()
}
