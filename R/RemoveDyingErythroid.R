#' RemoveDyingErythroid Function removes dying cells and erythroid cells
#'
#' @param ExprMat Expression matrix of genes X cells
#' @param Signature Signature genes to define a cell type for removal (default NULL)
#' @param MinPct.MT Minimal percentage of mitochondrial UMIs for a cell to be dying. default is 20%.
#' @param MinPct.Erythroid Minimal percentage of erythroid UMIs for a cell to be an erythrocyte. default is 20%.
#' @param MinPct.Signature Minimal percentage of signature UMIs for a cell to be a specific cell. default is 0.2%.
#'
#' @return an expression matrix after removing unwanted cells
#'
#' @importFrom Matrix readMM
#'
#' @export
#' @examples
#' RemoveDyingErythroid(data)
#'

RemoveDyingErythroid <- function(
  ExprMat,
  Signature = NULL,
  MinPct.MT = 0.2,
  MinPct.Erythroid = 0.2,
  MinPct.Signature = 0.002
){
  erythroid = c("HBG2","HBB","TFRC","TPT1","AHSP","FTL","HBG1","PRDX2","BLVRB","HBA1","HBA2","SLC4A1","ALAS2","HMGB2","HBD","FTH1","HBM","SLC25A37","UBB","SLC25A39","UCP2","SLC2A1","BSG","HMBS","HEMGN","TMCC2","RBM38","NCOA4","OAZ1","BCL2L1","FKBP8","PIM1","DCAF12","BPGM")
  geneset = match(erythroid, rownames(data))
  nErythroid <- colSums(as.matrix(ExprMat[geneset,]))

  #mitochondrial = c("MT-ND1","MT-ND2","MT-CO1","MT-CO2","MT-ATP8","MT-ATP6","MT-CO3","MT-ND3","MT-ND4L","MT-ND4","MT-ND5","MT-ND6","MT-CYB")
  geneset <- which(substr(rownames(data),1,3)=="MT-")
  if (length(geneset) < 2) {
    mitochondrial = c("ATP6","ATP8","COX1","COX2","COX3","CYTB","ND1","ND2","ND3","ND4","ND4L","ND5","ND6","RNR1","RNR2","TRNA","TRNC","TRND","TRNE","TRNF","TRNG","TRNH","TRNI","TRNK","TRNL1","TRNL2","TRNM","TRNN","TRNP","TRNQ","TRNR","TRNS1","TRNS2","TRNT","TRNV","TRNW","TRNY","MIR12136")
    geneset = match(erythroid, rownames(data))
  }
  nMT <- colSums(as.matrix(ExprMat[geneset,]))

  geneset = match(Signature, rownames(data))
  nSignature <- apply(ExprMat[geneset,], 2, sum)

  #megakaryocyte = c("ITGA2B","GP1BA","ITGB3","PF4","GP9","VWF","CXCL5","PPBP","TSPAN9","GP6")

  nUMIs = colSums(as.matrix(ExprMat))   # ==apply(ExprMat, 2, sum). Both don't work if too many cells/columns
  cells_selected = which(nMT/nUMIs < MinPct.MT & nErythroid/nUMIs < MinPct.Erythroid & nSignature/nUMIs < MinPct.Signature)
  return(ExprMat[, cells_selected])
}
