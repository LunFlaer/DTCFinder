#' Chromosomal positions of genes.
#'
#' A dataset containing genomic positons of over 20,000 genes
#'
#' @docType data
#' @usage data(gene_pos)
#' @format A data frame with over 20000 rows and 5 variables:
#' @keywords datasets
#' \describe{
#'   \item{ensembl_gene_id}{Ensembl gene ID, e.g. ENSG00000274847}
#'   \item{hgnc_symbol}{HGNC gene name, e.g. MAFIP}
#'   \item{chromosome_name}{chromosome, X --> 23, Y --> 24, MT --> 25}
#'   \item{start_position}{start position}
#'   \item{end_position}{end position}
#' }
#' @references  
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/nnnnnnnn}{PubMed})
#' @source \url{http://apr2020.archive.ensembl.org/}
#' @examples
#' data(gene_pos)
"gene_pos"
