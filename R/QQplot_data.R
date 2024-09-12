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
  Observed <- stats::rnorm(nrow(data))
  new_dat <- cbind(Observed, data[, -1])

  dat_plot <- new_dat %>%
    tidyr::gather(variable, value, -Observed)%>%
    dplyr::mutate(
      position = as.numeric(factor(variable, names(data)[-1])),
      order_col = (position - 1) %% 3
    ) %>%
    dplyr::group_by(order_col, position) |>
    dplyr::mutate(order_row = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::arrange(order_row, order_col) %>%
    dplyr::mutate(variable = factor(variable, levels = unique(variable)))

  dat_plot %>%
    ggplot2::ggplot(ggplot2::aes(sample = value, color = variable)) +
    ggplot2::stat_qq_line(
      col = "red",
      lwd = 0.5
    ) +
    ggplot2::theme(text = ggplot2::element_text(size = 14), legend.position = "none") +
    ggplot2::stat_qq() +
    ggplot2::facet_wrap(~variable, nrow=3) +
    ggplot2::ylab("Observed values") +
    ggplot2::xlab("Expected under normality")

  plotly::ggplotly()
}

