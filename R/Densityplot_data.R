#' Creating Density plot for a dataset
#'
#' @description Computes and draws kernel density estimate, which is a smoothed version of the histogram. 
#' This is a useful alternative to the histogram for continuous data that comes from an underlying smooth distribution. 
#' The `ggplot2` package is used here for creating the boxplot.
#'
#' @param data Proteomics expression dataset (original or normalized dataset) along with the protein information
#'
#' @return Interactive column-wise density plot
#'
#' @details This can also be used for comparing the original dataset with the normalized dataset.
#'
#' @seealso
#' `geom_density()`
#'
#' @export
#'
#' @examples Densityplot_data(yeast_data)
#' @examples Densityplot_data(rlr_knn_yeast_data)
Densityplot_data <- function (data){

  new_data <- as.data.frame(data)
  new_data <- new_data[,-1]
  
  dat_plot <- new_data %>%
    tidyr::gather(variable, value) %>%
    dplyr::mutate(
      position = as.numeric(factor(variable, names(data)[-1])),
      order_col = (position - 1) %% 3
    ) %>%
    dplyr::group_by(order_col, position) |>
    dplyr::mutate(order_row = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::arrange(order_row, order_col) %>%
    dplyr::mutate(variable = factor(variable, levels = unique(variable)))

  density_plot <- dat_plot %>%
    ggplot2::ggplot(ggplot2::aes(x=value) ) +
    ggplot2::geom_density (fill= "#69b3a2") +
    ggplot2::facet_wrap(~variable, scales="free", nrow = 3)+
    ggplot2::theme_gray()+
    ggplot2::theme(text = ggplot2::element_text(size = 14))
  suppressWarnings(return(plotly::ggplotly(density_plot)))
}
