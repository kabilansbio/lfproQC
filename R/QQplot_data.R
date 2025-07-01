#' Creating QQ-Plot for a dataset
#'
#' @description A Q–Q plot (quantile-quantile plot) is a plot of the quantiles of two distributions against each other, or a plot based on estimates of the quantiles. 
#' The normality of the data can be understand by this plot.
#'
#' @param data Proteomics expression dataset (original or normalized dataset)
#'
#' @details This can be used for comparing the original dataset with the
#' normalized dataset.
#'
#' @return Interactive column-wise QQ-plot
#' @export
#'
#' @examples qqplot <- QQplot_data(rlr_knn_yeast_data)
QQplot_data <- function(data) {
  Observed <- feature <- ppoints <- theoretical <- NULL
  data <- as.data.frame(data)
  colnames(data) <- make.unique(colnames(data))
  
  new_dat <- data[, -1, drop = FALSE]
  new_dat$Observed <- stats::rnorm(nrow(new_dat))
  
  dat_plot <- tidyr::pivot_longer(new_dat, cols = -Observed, names_to = "feature", values_to = "value")
  
  dat_plot <- dat_plot %>%
    dplyr::group_by(feature) %>%
    dplyr::mutate(
      theoretical = stats::qnorm(ppoints(length(value))),
      value = sort(value)
    ) %>%
    dplyr::ungroup()
  
  p <- ggplot2::ggplot(dat_plot, ggplot2::aes(x = theoretical, y = value, color = feature)) +
    ggplot2::geom_point(size = 1, alpha = 0.6) +
    ggplot2::geom_abline(slope = 1, intercept = 0, col = "red", lwd = 0.5) +
    ggplot2::facet_wrap(~feature, nrow = 3, scales = "free") +
    ggplot2::ylab("Observed values") +
    ggplot2::xlab("Theoretical quantiles") +
    ggplot2::theme(text = ggplot2::element_text(size = 14), legend.position = "none")
  
  plotly::ggplotly(p)
}

