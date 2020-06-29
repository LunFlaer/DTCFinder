#' CNV_Simulator Function infers CNV profiles of cells in terms of posterio probability  
#' based on a fitted two-component Gaussian Mixture Model over the log2(CPM/10+1) expression matrix
#'
#' @param ExprMat A genes X cells expression matrix of log2(CPM/10+1) scaled scRNA-seq counts
#' @param minUMIs Minimal UMI counts for a gene to be kept over all cells (default = 10)
#' @param minCells Minimal expressed cell number for a gene to be kept over all cells (default = 1)
#' @param minGenes Minimal region size in terms of genes number within the region (default = 100)
#' @param repetitions Loop number of fitting to a GMM to obtain an optimal Gaussian mixture model
#'
#' @return a matrix of posterior probabilities as CNV estimates
#'
#' @importFrom edgeR cpm
#' @importFrom mixtools normalmixEM
#'
#' @export
#' @examples
#' CNV_Simulator(data)
#'

CNV_Simulator <- function(
  ExprMat,
  minUMIs = 10,
  minCells = 1,
  minGenes = 100,
  repetitions = 8
){
  data = log2(edgeR::cpm(ExprMat, lib.size = NULL, log = F, prior.count = 0) / 10 + 1)  #	Counts per million for ExprMat, scaled by 10, and transformed with pseudo-count = 1 and log2
  data(gene_pos)
  data = data[intersect(rownames(data), gene_pos[,"hgnc_symbol"]), ]
  
  
  data = data[which(rowSums(data) > minUMIs), ]
  nCells = apply(data, 1, function(x){sum(x>0, na.rm=TRUE)})
  data = data[nCels > minCells, ]
  normFactor = colMeans(data)

  gene_positions = gene_pos[sort(match(rownames(data), gene_pos[,2])),]
  gene_positions = gene_positions[!(is.na(suppressWarnings(as.integer(gene_positions[,3])))),]
  gene_positions = gene_positions[which(as.integer(gene_positions[,3])>0), ]
  gene_positions = gene_positions[which(as.integer(gene_positions[,3])<25), ]

  regions = c()
  n = 1
  i = 1
  chrom = as.integer(gene_positions[i,3])
  start = as.integer(gene_positions[i,4]) - 1

  while (i < (nrow(gene_positions) - RegionSize)) {
    for (j in 1:RegionSize) {
      if (is.na(as.integer(gene_positions[i+j,3]))) {
        chr = 0
      }
      else {
        chr = as.integer(gene_positions[i+j,3])
      }
      if (chrom !=  chr) {
        break
      }
    }
    if (j < RegionSize) {
      j = j-1
    }
    i = i+j+1
    end = as.integer(gene_positions[i-1,5]) + 1
    reg_new = c(chrom = chrom, start = start, end = end, length = (end - start))
    regions = rbind(regions, reg_new)
    rownames(regions)[n] <- n
    n = n+1
    if (i<nrow(gene_positions)) {
      if (is.na(as.integer(gene_positions[i,3]))) {
        chrom = 0
        break
      }
      else {
        chrom = as.integer(gene_positions[i,3])
      }
      start = as.integer(gene_positions[i,4]) - 1
    }
  }

  GenePositions <- data.frame(ensembl_gene_id=gene_positions[,1], hgnc_symbol=gene_positions[,2], chromosome_name=as.integer(gene_positions[,3]), start_position=as.integer(gene_positions[,4]), end_position=as.integer(gene_positions[,5]))
  MMPP = c()
  for(i in 1:nrow(regions)){
    genes_region = GenePositions[which(GenePositions[,3]==regions[i,1] & GenePositions[,4]>regions[i,2] & GenePositions[,5]<regions[i,3]), 2]
    if(length(genes_region) >= minGenes){
      expr_region = scale(colMeans(data[intersect(genes_region, row.names(data)),]) - normFactor)
	    bestlog = (-Inf)
	    bestGMM = NULL
	    for (i in 1:repetitions){
		    print(paste0("Fitting GMM for region ", chr, ":", regions[i,2], "-", regions[i,3]," iteration ",i))
		    GMM = tryCatch(mixtools::normalmixEM(expr_region, k = 2, maxit = 1000, maxrestarts = 10), 
          error = function(e) {print(paste0("EM algorithm did not converge for region ", chr, ":", regions[i,2], "-", regions[i,3])); GMM = NULL})
		    if(!is.null(GMM)){
			    if (GMM$loglik > bestlog){
			      bestlog = GMM$loglik
				    bestGMM = GMM
			    }
		    }
	    }
      if(!is.null(bestGMM)){
        if (bestGMM$mu[1] > bestGMM$mu[2]){
          pp = bestGMM$posterior[,1]
        }
        else{
          pp = bestGMM$posterior[,2]
        }
        MMPP = cbind(MMPP, pp)
        colnames(MMPP)[ncol(MMPP)] = rownames(regions)[i]
      }
    }
  }
  return(MMPP)
}
