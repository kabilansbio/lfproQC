#' Creating Correlation matrix plot for a dataset
#'
#' @description A graphical display of a correlation matrix.
#'
#' @param data Proteomics expression dataset (original or normalized dataset) along with the protein information
#'
#' @details This can also be used for comparing the original dataset with the normalized dataset.
#'
#' @return Interactive corrleation matrix plot
#' @export
#'
#' @examples Corrplot_data(yeast_data)
#' @examples Corrplot_data(rlr_knn_yeast_data)
Corrplot_data <- function(data){
  x <- data[,-1]
  y <- as.matrix(x)
  rt <- Hmisc::rcorr(y)
  mtlr <- suppressWarnings(reshape::melt(rt$r))
  mtlp <- suppressWarnings(reshape::melt(rt$P))
  p.value <- mtlp$value
  gx <- ggplot2::ggplot(mtlr, ggplot2::aes(X1, X2, fill = value, label=p.value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "cyan",  high = "red")
  plotly::ggplotly(gx)
}

