#' Creating MDS plot for a dataset
#'
#' @description Multi-dimensional scaling (MDS) plots showing a 2-dimensional projection of distances between the dataset samples.
#'
#' @param data Normalized and imputed Proteomics expression dataset along with the protein information
#'
#' @return MDS plot
#'
#' @seealso
#' `mdsPlot`
#'
#' @export
#'
#' @examples MDSplot_data(rlr_knn_yeast_data)
MDSplot_data <- function (data) {
  X <- Y <- Sample <- NULL
  data1 <- data[,-1]
  num_data <- data.frame(sapply(data1, function(x) as.numeric(as.character(x))))
  log2.data.matrix <- log2(num_data)
  log2.distance.matrix <- matrix(0,
                                 nrow=ncol(log2.data.matrix),
                                 ncol=ncol(log2.data.matrix),
                                 dimnames=list(colnames(log2.data.matrix),
                                               colnames(log2.data.matrix)))
  for(i in 1:ncol(log2.distance.matrix)) {
    for(j in 1:i) {
      log2.distance.matrix[i, j] <-
        mean(abs(log2.data.matrix[,i] - log2.data.matrix[,j]))
    }
  }

  mds.stuff <- stats::cmdscale(stats::as.dist(log2.distance.matrix),
                        eig=TRUE,
                        x.ret=TRUE)

  mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)

  mds.values <- mds.stuff$points
  mds.data <- data.frame(Sample=rownames(mds.values),
                         X=mds.values[,1],
                         Y=mds.values[,2])

  mds_plot <- ggplot2::ggplot(data=mds.data, ggplot2::aes(x=X, y=Y, label=Sample)) +
    ggplot2::geom_text() +
    ggplot2::theme_bw() +
    ggplot2::xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
    ggplot2::ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
    ggplot2::ggtitle("MDS plot using avg(logFC) as the distance")+
    ggplot2::theme(text = ggplot2::element_text(size=14))+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  return(plotly::ggplotly(mds_plot))
}
