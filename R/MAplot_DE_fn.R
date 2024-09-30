#' Find out the Up and Down regulated proteins from MA plot
#'
#' @description MA plot is used for visualizing the differentially expressed proteins by plotting the log mean intensity data in x axis and log fold change values in y axis.
#'
#' This function can be used for visualizing the up regulated, down regulated, and non-significant proteins along with their information.
#'
#' @param top_table Top table information
#' @param x1 Cut-off limit for down-regulated proteins
#' @param x2 Cut-off limit for up-regulated proteins
#' @param p Cut-off limit for p-values
#'
#' @return
#' `Result` Top table along with up, down, significant and non-significant protein information.
#'
#' `MA plot` Interacive MA plot with the details of up and down regulated proteins
#'
#' `Up-regulated` Up-regulated protein information
#'
#' `Down-regulated` Down-regulated protein information
#'
#' `Non-significant` Non-significant protein information
#'
#' @export
#'
#' @examples result <- MAplot_DE_fn(yeast_top_table, -1, 1, 0.05)
#' @examples result$`MA Plot`
#' @examples result$`Result`
#' @examples result$`Up-regulated`
#' @examples result$`Down-regulated`
#' @examples result$`Non-significant`
MAplot_DE_fn <- function(top_table, x1=NULL, x2=NULL, p=NULL){

  AveExpr <- logFC <- diffexpressed <- NULL

  if(is.null(x1) && is.null(x2) && is.null(p)){
    AveExpr <- logFC <- NULL
    ID <- rownames(top_table)
    plot1 <- plotly::ggplotly(ggplot2::ggplot (top_table, ggplot2::aes(x = AveExpr, y = logFC, label = ID))+
      ggplot2::geom_point())
    return(plot1)
  }

  else{
  #Changing the first column as a row name
  res_table <- top_table

  #Add a new variable for three level factors - Up, Down and NO

  res_table$diffexpressed <- "Non-significant"

  res_table$diffexpressed[res_table$logFC >=x2 & res_table$adj.P.Val <= p] <- "Up"

  res_table$diffexpressed[res_table$logFC <=x1 & res_table$adj.P.Val <= p] <- "Down"

  #res_table$diffexpressed[res_table$logFC>x2 & res_table$logFC<x1 & res_table$adj.P.Val < p] <- "Significant"

  res_table$diffexpressed[res_table$logFC<x1 & res_table$logFC>x2 & res_table$adj.P.Val >p] <- "Non-signficant"


  #Create the plot
  ID <- rownames(res_table)
  plot <- ggplot2::ggplot(data = res_table, ggplot2::aes(x=AveExpr, y=logFC,
                                                         label = ID, col = diffexpressed)) +
    ggplot2::geom_point()+
    ggplot2::theme_minimal()+
    ggplot2::scale_color_manual(values=c("Non-significant" = "black", "Down" = "blue", "Up" = "green"))+
    ggplot2::geom_hline(yintercept = c(x1,x2), col = "azure4", linewidth = 1)+

    ggplot2::ggtitle("MA Plot")+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  #Up-regulated proteins
  Up_protein <- subset(res_table, diffexpressed == "Up")

  #Down-regulated proteins
  Down_protein <- subset (res_table, diffexpressed == "Down")

  #Non-significant proteins
  Nonsignificant_protein <- subset (res_table, diffexpressed == "Non-significant")

  result = list ("Result" = res_table, "MA Plot" =plotly::ggplotly(plot),
                 "Up-regulated" = Up_protein, "Down-regulated" = Down_protein,
                 "Non-significant" = Nonsignificant_protein)
  return(result)
  }
}

