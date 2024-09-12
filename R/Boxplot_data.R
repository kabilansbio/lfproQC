#' Creating Boxplot for a dataset
#'
#' @description The box and whiskers plot displays the distribution of a continuous variable. 
#' It visualises five summary statistics (the median, two hinges and two whiskers), and all "outlying" points individually. 
#' The `ggplot2` package is used here for creating the boxplot.
#'
#' @param data Proteomics expression dataset (original or normalized dataset)
#'
#' @return Interactive box and whiskers plot
#'
#' @details This can also be used for comparing the original dataset with the normalized dataset.
#'
#' @seealso
#' `geom_boxplot()`
#'
#' @export
#'
#' @examples Boxplot_data(yeast_data)
#' @examples Boxplot_data(rlr_knn_yeast_data)
Boxplot_data <- function (data){
  variable <- value <- NULL

  meltData <- reshape2::melt(data)
  plot <- ggplot2::ggplot(meltData, ggplot2::aes(x = variable, y = value))+
    ggplot2::geom_boxplot(ggplot2::aes(fill=variable),outlier.size = 0.8, show.legend = FALSE)+
    ggplot2::labs(x = "Columns")+
    ggplot2::ggtitle("Boxplot")+
    ggplot2::theme(text = ggplot2::element_text(size=16))+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  return(plotly::ggplotly(plot))
}
