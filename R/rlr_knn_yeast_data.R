#' Normalized and imputed complete yeast lysate - UPS1 benchmark dataset
#'
#' This is the groupwise normalized and missing values imputed dataset of the complete yeast lysate - UPS1 benchmark dataset. 
#' Normalization has been done by RLR normalization method and missing values imputation has been done by KNN imputation method.
#'
#' @format A data frame with 835 rows and 7 variables:
#' \describe{
#'   \item{Majority protein IDs}{Protein ID information}
#'   \item{A1}{1st sample group, 1st technical replicate}
#'   \item{A2}{1st sample group, 2nd technical replicate}
#'   \item{A3}{1st sample group, 3rd technical replicate}
#'   \item{B1}{2nd sample group, 1st technical replicate}
#'   \item{B2}{2nd sample group, 2nd technical replicate}
#'   \item{B3}{2nd sample group, 3rd technical replicate}
#'}
"knn_rlr_yeast_data"
