#' Creating the top table
#'
#' @description Top table can be used for identifying the pairwise differential abundance analysis of proteins in the dataset.
#'
#' @param data Normalized and missing values imputed expression dataset containing protein information
#' @param groups Group information about the input data
#' @param ch_gr1 Group number of the dataset for pairwise comparison with the another group
#' @param ch_gr2 Group number of the dataset to be compared with the chosen group
#'
#' @return Top table consists of following values
#'
#' `logFC` -  Log fold change values,
#'
#' `AveExpr` - Average intensity values,
#'
#' `t`- t-statistic values,
#'
#' `P.Value` - P-values,
#'
#' `adj.P.Val` - Adjusted P-values,
#'
#' `B`- B-statistic values
#'
#' @seealso `limma::topTable`
#'
#' @export
#'@examples
#' \donttest{
#' top_table <- top_table_fn(rlr_knn_yeast_data, yeast_groups, 2, 1)
#' top_table
#' }
top_table_fn <- function(data, groups, ch_gr1, ch_gr2){
  #Changing the first column as a row name
  data2 <- as.data.frame(data[,-1])
  Protein_info <- as.data.frame(data[,1])
  colnames(Protein_info)[1]<- "Protein_det"
  exp.data1 <- cbind(Protein_info, data2)
  rownames(exp.data1) <- make.names(exp.data1$Protein_det, unique = TRUE)
  exp.data <- exp.data1[,-1]


  #Creating the model matrix
  eset <- as.matrix(exp.data)
  Groups <- groups$Groups
  group_name <- unique(Groups)
  cont.grp <- paste(group_name[ch_gr1], "-", group_name[ch_gr2], sep = "")
  design <- stats::model.matrix(~0+factor(Groups))
  colnames(design)<- group_name

  contr.matrix <- limma::makeContrasts(contrasts = cont.grp, levels=design)

  #Fitting the model
  fit <-limma::lmFit(eset, design)
  fit <- limma::contrasts.fit(fit, contrasts = contr.matrix)
  fit <- limma::eBayes(fit)
  res <- limma::topTable(fit, number= Inf, adjust.method = "BH")
  return(res)
}

