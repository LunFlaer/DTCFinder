#' MalignantCellsFinder Function identifies putative malignant cells based on inferred CNV profiles of cells
#' using two-component Gaussian Mixture Models fitted to 10X single cells expression profiles
#'
#' @param PathToData Folder with barcodes.tsv, features.tsv (or genes.tsv), matrix.mtx
#' @param pvalue P value to define outliers relative to CD45-positive cells. Default = 0.05
#' @param p_adj P value to define immune cells and malignant cells. Default = 1e-6
#'
#' @return a vector of predicted malignant cells
#'
#' @importFrom Matrix readMM
#' @importFrom edgeR cpm
#' @importFrom mixtools normalmixEM
#' @importFrom ggplot2 ggplot
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#' @examples
#' MalignantCellsFinder("outs/filtered_feature_bc_matrix/")
#'

MalignantCellsFinder <- function(
  PathToData = ".",
  pvalue = 0.05,
  p_adj = 1e-6  # 1e-2 for tumor tissue, such as case O
){
  ExprMat <- FilterData(PathToData)
  ExprMat <- RemoveDyingErythroid(ExprMat, Signature = NULL, MinPct.MT = 0.2, MinPct.Erythroid = 0.2, MinPct.Signature = 0.002)
  CD45 = ExprMat["PTPRC",]
  CD45p_cells = which(CD45>0)
  CD45n_cells = which(CD45==0)
  if (length(CD45p_cells) < 1) {
      stop("Error: No CD45-positive cells found. Exit.")
  }

  MMPP <- CNV_Simulator(ExprMat)

  distribution = data.frame(average=apply(MMPP[CD45p_cells,],2,mean),stdev=apply(MMPP[CD45p_cells,],2,sd))
  low = apply(distribution,1,function(x) qnorm(pvalue, mean = x[1], sd = x[2], lower.tail = TRUE, log.p = FALSE))
  up = apply(distribution,1,function(x) qnorm(pvalue, mean = x[1], sd = x[2], lower.tail = FALSE, log.p = FALSE))

  outlier_count = rep(0, nrow(MMPP))
  for (i in 1:nrow(MMPP)) {
    outlier_count[i] = sum(MMPP[i,]<low | MMPP[i,]>up)
  }

  his_CD45p = rep(0, ncol(MMPP)+1)
  his_CD45n = rep(0, ncol(MMPP)+1)
  for (i in 1:(ncol(MMPP)+1)) {
    his_CD45p[i] = sum(outlier_count[CD45p_cells]==(i-1))
    his_CD45n[i] = sum(outlier_count[CD45n_cells]==(i-1))
  }

  mu = mean(outlier_count[CD45p_cells])
  sigma = sd(outlier_count[CD45p_cells])
  Amplitute = mean(head(sort(his_CD45p,decreasing=TRUE),3))

  names(his_CD45p) <- paste(seq_along(his_CD45p)-1)
  bar <- data.frame(x = seq_along(his_CD45p)-1, y = his_CD45p)
  model <- nls(y ~ A*exp(-0.5*((x-xc)/w)^2), start = c(xc=mu, w=sigma, A=Amplitute) , data = bar)
  r_square <- RSquared(model)
  if (r_square[2] < 0.8) {
      stop(paste0("Error: Exceptional Regions of CD45-positive cells do not follow a Guasian distribution well: adjusted R-square = ", r_square[2]))
  }

  v <- summary(model)$parameters[,"Estimate"]

  prob = rep(0, ncol(MMPP)+1)
  for (n in 0:ncol(MMPP)) {
    prob[n+1] = pnorm(n, mean = v[1], sd = v[2], lower.tail = FALSE, log.p = FALSE)
  }

  #library("ggplot2")
  png(filename = "ExceptionalRegionDistribution.png", width = 3300, height = 1800, pointsize = 12, res=600)
  ggplot2::ggplot(data=bar, aes(x=x, y=y)) + xlab("Number of Outlier Regions") + ylab("Cell Count") +
    geom_bar(stat="identity", fill="steelblue") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_function(fun = function(x) dnorm(x, mean = v[1], sd = v[2]) * v[3] * 2*v[2]*sqrt(pi/2), color = "darkred", size = 0.5) +
    ggtitle("Exceptional Regions Distribution of CD45-positive Cells")
  dev.off()

  CD45p_cells_adj = setdiff(CD45p_cells,which(outlier_count>=(head(which(prob<p_adj),1)-1)))  
  distribution_adj = data.frame(average=apply(MMPP[CD45p_cells_adj,],2,mean),stdev=apply(MMPP[CD45p_cells_adj,],2,sd))
  low = apply(distribution_adj,1,function(x) qnorm(p, mean = x[1], sd = x[2], lower.tail = TRUE, log.p = FALSE))
  up = apply(distribution_adj,1,function(x) qnorm(p, mean = x[1], sd = x[2], lower.tail = FALSE, log.p = FALSE))

  outlier_count_adj = rep(0, nrow(MMPP))
  for (i in 1:nrow(MMPP)) {
    outlier_count_adj[i] = sum(MMPP[i,]<low | MMPP[i,]>up)
  }

  his_CD45p_adj = rep(0, ncol(MMPP)+1)
  his_CD45n_adj = rep(0, ncol(MMPP)+1)
  his_CD45n = rep(0, ncol(MMPP)+1)
  for (i in 1:(ncol(MMPP)+1)) {
    his_CD45p_adj[i] = sum(outlier_count_adj[CD45p_cells_adj]==(i-1))
    his_CD45n[i] = sum(outlier_count_adj[CD45n_cells]==(i-1))
  }

  mu_adj = mean(outlier_count_adj[CD45p_cells_adj])
  sigma_adj = sd(outlier_count_adj[CD45p_cells_adj])
  Amplitute_adj = mean(head(sort(his_CD45p_adj,decreasing=TRUE),3))

  names(his_CD45p_adj) <- paste(seq_along(his_CD45p_adj)-1)
  bar_adj <- data.frame(x = seq_along(his_CD45p_adj)-1, y = his_CD45p_adj)
  model_adj <- nls(y ~ A*exp(-0.5*((x-xc)/w)^2), start=c(xc=mu_adj, w=sigma_adj, A=Amplitute_adj) , data = bar_adj)
  v <- summary(model_adj)$parameters[,"Estimate"]

  prob_adj = rep(0, ncol(MMPP)+1)
  for (n in 0:ncol(MMPP)) {
    prob_adj[n+1] = pnorm(n, mean = v[1], sd = v[2], lower.tail = FALSE, log.p = FALSE)
  }

  CD45p_cells_adj = setdiff(CD45p_cells_adj,which(outlier_count_adj>=(head(which(prob_adj<p_adj),1)-1)))
  CD45n_cells_adj = setdiff(CD45n_cells,which(outlier_count_adj<=(head(which(prob_adj<(p_adj/100)),1)-1)))
  for (i in 1:(ncol(MMPP)+1)) {
    his_CD45p_adj[i] = sum(outlier_count_adj[CD45p_cells_adj]==(i-1))
    his_CD45n_adj[i] = sum(outlier_count_adj[CD45n_cells_adj]==(i-1))
  }

  bin = length(his_CD45p_adj)
  types = rep("immune", n)
  types[his_CD45n_adj>0] <- "tumor-like"
  while (bin > 0) {
    if (his_CD45n_adj[bin]>0 | his_CD45p_adj[bin]>0) {
      break
    }
    bin = bin - 1
  }

  total_CD45p_cells = length(CD45p_cells_adj)
  total_CD45p_cells = length(CD45n_cells_adj)
  if (total_CD45p_cells == 0) total_CD45p_cells = 1
  bar_adj <- data.frame(x = (seq_along(his_CD45p_adj)-1)[1:bin], yp = his_CD45p_adj[1:bin]/total_CD45p_cells, yn = his_CD45n_adj[1:bin]/total_CD45n_cells, celltypes = types[1:bin])

  colors = c(rev(RColorBrewer::brewer.pal(9, "YlGnBu"))[1:9], RColorBrewer::brewer.pal(9, "YlOrRd")[2:9])
  png(filename = "MalignantCellsDistribution.png", width = 3300, height = 1800, pointsize = 9, res=600)
  ggplot(data=bar_adj, aes(x=x, fill=celltypes)) + xlab("Number of Outlier Regions") + ylab("Relative Cell Count") +
  geom_bar(aes(y=yp), stat="identity", show.legend = TRUE) +
  geom_bar(aes(y=yn), stat="identity", show.legend = TRUE) +
  scale_fill_manual(values=c(colors[4], colors[13])) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_function(fun = function(x) dnorm(x, mean = v[1], sd = v[2]) * v[3] * 2*v[2]*sqrt(pi/2), color = colors[16], size = 0.5) +
  theme(text = element_text(size=12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position="bottom") +
  ggtitle("Exceptional Regions Distribution")
  dev.off()

  return(CD45n_cells_adj)
}
