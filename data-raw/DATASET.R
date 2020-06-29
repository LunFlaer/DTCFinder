## code to prepare `gene_pos` dataset
#' 
#' getGenePositions() Receive genomic coordinates of a gene list
#'
#' This function retrieves the Ensembl database for chromosomal positions of HUGO genes
#' @param GeneNames A vector of gene names in HUGO format.
#' @param EnsemblVersion Version of the ENSEMBL database. Defaul: apr2020 v100 (GRCh38p13).
#' @param Species Either human or mouse
#' @keywords Chromosomal positions
#' @importFrom biomaRt useMart getBM
#' @export
#' @examples
#' getGenePositions(GeneNames=c("EGFR","PDGFRA"))
#'
getGenePositions = function(GeneNames, EnsemblVersion="apr2020.archive.ensembl.org", Species="human"){
  if (Species=="human"){
    ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = EnsemblVersion)
    gene_positions <- biomaRt::getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'), filters = 'hgnc_symbol', values = GeneNames, mart = ensembl)
  }
  else {
    ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = EnsemblVersion)
    gene_positions <- biomaRt::getBM(attributes=c('ensembl_gene_id','mgi_symbol','chromosome_name','start_position','end_position'), filters = 'mgi_symbol', values = GeneNames, mart = ensembl)
  }
  gene_positions=gene_positions[!duplicated(gene_positions[,2]),]
  gene_positions[which(gene_positions[,3]=="X"),3] = 23
  gene_positions[which(gene_positions[,3]=="Y"),3] = 24
  gene_positions[which(gene_positions[,3]=="MT"),3] = 25
  gene_positions[which(nchar(gene_positions[,3]) > 2),3] = 0
  gene_positions = gene_positions[order(as.numeric(gene_positions[,3])*1e9 + as.numeric(gene_positions[,4]),decreasing=F),]
  return(gene_positions)
}

gene_pos = getGenePositions(rownames(data))
usethis::use_data(gene_pos, overwrite = TRUE)

